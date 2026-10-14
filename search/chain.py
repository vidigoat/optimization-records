"""A search that works in turns on one line of attack, instead of many one-shot programs.

The batch search (common.evolve) asks for many independent programs and keeps the best. That is
wide but shallow: no program ever sees the output of its own experiment. Here one chain of model
calls works on one direction. Each turn the model writes a program, we run it, and the next turn
it reads what the program printed, the verifier's verdict, and the notes it left for itself. A
program may be an experiment and need not be a final answer.

Still one model call at a time per chain, each reserved in the ledger before it is sent.
"""
import json
import os
import re
import threading
import time

from . import budget, llm, sandbox

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

RULES = """## How this lab works
You work in turns. Each turn you write one Python program. It runs for at most {seconds} seconds on
{cores} with no network, and on the next turn you are shown the last 3000 characters it printed to
stdout, the last 1500 of stderr, and the official verdict on its answer. You remember nothing
between turns except the NOTES you write, so put there what you learned, what is ruled out and
what to do next. A turn may be an experiment: print measurements that will tell you what to try.
You may use C for speed: write the source to a file in the current directory, compile it with
`cc -O3 -march=native`, and run it with subprocess. Standard library only otherwise. {parallel} Print the answer as one line of JSON on the last line of stdout, and print the best
answer found so far before the time limit, because a program stopped at the limit reports nothing more.
"""

REPLY = """## Your turn ({turn} of {turns})
Reply with one line beginning `IDEA:` saying what this turn tests, then one ```python block, then
a line `NOTES:` followed by your updated notes (at most 400 words; they replace the old notes).
Keep your reasoning short so the reply is not cut off before the program and the notes."""


def _state(path):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return {"turns": [], "notes": "", "best": None, "usd": 0.0}


def _prompt(task, st, angle, turns):
    cores = getattr(task, "cores", 1)
    parts = [task.brief, RULES.format(
        seconds=task.run_seconds, cores="one core" if cores == 1 else f"{cores} cores",
        parallel="Use a single process." if cores == 1 else
        f"Use all {cores} cores: start {cores} copies of your compiled search with different seeds or different parts of the "
        "search space, and merge their results.")]
    if angle:
        parts.append(f"## Your direction\n{angle}\nStay with this direction and take it as deep as it goes.\n")
    if st["notes"]:
        parts.append(f"## Your notes from earlier turns\n{st['notes']}\n")
    if st["best"]:
        b = st["best"]
        parts.append(f"## Best program of this chain so far: score {b['score']} ({b['summary']})\n```python\n{b['program']}\n```\n")
    if st["turns"]:
        last = st["turns"][-1]
        if not st["best"] or last["program"] != st["best"]["program"]:
            parts.append(f"## Last turn's program\n```python\n{last['program']}\n```\n")
        parts.append(f"## What last turn's program printed\nstdout (tail):\n{last['stdout']}\nstderr (tail):\n{last['stderr']}\n"
                     f"Verdict: score {last['score']} - {last['summary']}\n")
        parts.append("## All turns so far\n" + "\n".join(
            f"- turn {t['turn']}: score {t['score']} ({t['summary'][:110]}) - {t['idea'][:150]}" for t in st["turns"]) + "\n")
    parts.append(REPLY.format(turn=len(st["turns"]) + 1, turns=turns))
    return "\n".join(parts)


def run_chain(task, name, turns, *, angle="", effort="high", max_tokens=32000, max_usd=10.0, slots=None,
              stop_at=None, model=llm.DEFAULT_MODEL):
    """Run one chain until `turns` turns, its own spending limit, the goal's cap, or `stop_at`."""
    path = os.path.join(ROOT, "runs", f"chain_{name}.json")
    st = _state(path)
    empty = 0
    while len(st["turns"]) < turns and st["usd"] < max_usd:
        n = len(st["turns"]) + 1
        try:
            text, usd = llm.ask(task.name, _prompt(task, st, angle, turns), model=model, effort=effort,
                                max_tokens=max_tokens, note=f"chain {name} turn {n}")
        except budget.OverBudget as e:
            print(f"[{name}] stopping: {e}", flush=True)
            break
        except Exception as e:
            print(f"[{name}] turn {n}: API error {type(e).__name__}: {str(e)[:120]}", flush=True)
            time.sleep(60)
            continue
        st["usd"] += usd
        if not text.strip() and empty < 2:
            # The API sometimes declines a harmless prompt outright. Ask again before spending a turn on it.
            empty += 1
            print(f"[{name}] turn {n}: empty reply, asking again", flush=True)
            continue
        empty = 0
        idea = (re.search(r"IDEA:\s*(.+)", text) or [None, "(no idea line)"])[1][:300]
        blocks = re.findall(r"```python\n(.*?)```", text, re.S)
        program = max(blocks, key=len) if blocks else ""
        notes = text.rsplit("NOTES:", 1)[1].strip()[:4000] if "NOTES:" in text else ""
        if notes:
            st["notes"] = notes
        if program:
            with (slots or threading.Lock()):
                answer, out, err = sandbox.run_full(program, task.run_seconds)
            score, summary = task.score(answer)
        else:
            answer, out, err, score, summary = None, "", "", 0, "no program in reply"
        st["turns"].append({"turn": n, "idea": idea, "program": program, "stdout": out, "stderr": err,
                            "score": score, "summary": summary, "usd": round(usd, 4), "t": time.time()})
        if program and score > 0:
            with open(os.path.join(ROOT, "runs", f"chain_{name}_answer_{n:03d}_score{score}.json"), "w") as f:
                json.dump(answer, f)
            if not st["best"] or score > st["best"]["score"]:
                st["best"] = {"score": score, "summary": summary, "program": program}
        tmp = path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(st, f, indent=1)
        os.replace(tmp, path)
        print(f"[{name}] turn {n}: score {score} ({summary[:90]})  ${usd:.2f} (chain ${st['usd']:.2f})  | {idea[:140]}", flush=True)
        if stop_at is not None and score >= stop_at:
            print(f"[{name}] TARGET {stop_at} REACHED in turn {n}", flush=True)
            break
    return st


def run_chains(task, angles, turns, *, prefix, parallel=2, **kw):
    """Several chains at once, one per line of attack. Programs share `parallel` CPU slots."""
    slots = threading.BoundedSemaphore(parallel)
    threads = [threading.Thread(target=run_chain, args=(task, f"{prefix}{i}", turns),
                                kwargs=dict(angle=a, slots=slots, **kw)) for i, a in enumerate(angles)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
