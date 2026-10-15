"""The spend ledger. Every API call goes through here, and it fails closed.

The whole project has a hard cap (750 USD) split into per-goal caps. Before a
call we reserve its worst-case cost; if that would cross the goal's cap or the
project cap, the call is refused. After the call we record what it really cost
from the usage the API reports. The ledger is an append-only file on disk, so
a crash or a second process cannot lose track of money already spent.
"""
import fcntl
import json
import os
import threading
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.environ.get("MOONSHOTS_LEDGER", os.path.join(ROOT, "runs", "ledger.jsonl"))

PROJECT_CAP = 750.00
# What the key in use can actually pay for. The credits sit in two accounts, so the first key
# is capped below the project total; raise this when the second key takes over.
ACCOUNT_CAP = float(os.environ.get("MOONSHOTS_ACCOUNT_CAP", "470"))
# Two keys from two accounts. Key 1 is used until its cap, then key 2. Key 2's account holds a $380.40
# credit and a separate $100 credit the owner wants left alone, so its cap stays below $380.
ACCOUNT_CAPS = {"1": ACCOUNT_CAP, "2": float(os.environ.get("MOONSHOTS_ACCOUNT2_CAP", "370"))}
# Shannon, Erdős and Goldbach have ended; their caps were lowered to what they spent, and the freed
# money split between the three record-table goals. The project total stays 750.
GOAL_CAPS = {"shannon": 200.00, "erdos": 28.00, "machines": 51.00, "reserve": 40.00,
             "packing": 20.00, "matmul": 70.00, "tao": 211.00,     # matmul ended without a result; $70 moved to tao
             "codes": 130.00}                                      # packing finished; its unused $130 moved to codes

# USD per million tokens: (input, output, cache read, 5-minute cache write).
# From platform.claude.com/docs/en/about-claude/pricing, fetched 13 Oct 2026.
PRICES = {
    "claude-fable-5-1":  (10.00, 50.00, 0.25, 12.50),
    "claude-opus-5-5":   (4.00, 20.00, 0.20, 5.00),
    "claude-sonnet-5-5": (2.00, 10.00, 0.10, 2.50),
    "claude-haiku-5-5":  (0.10, 0.50, 0.01, 0.125),   # prompts up to 100k tokens
}
BATCH_DISCOUNT = 0.5


class OverBudget(RuntimeError):
    pass


# Parallel rounds share one process: checking the cap and claiming the money must be one step.
_RESERVE_LOCK = threading.Lock()


def cost(model, input_tokens=0, output_tokens=0, cache_read=0, cache_write=0, batch=False):
    pin, pout, pread, pwrite = PRICES[model]
    usd = (input_tokens * pin + output_tokens * pout + cache_read * pread + cache_write * pwrite) / 1e6
    return usd * (BATCH_DISCOUNT if batch else 1.0)


def _entries():
    if not os.path.exists(LEDGER):
        return []
    with open(LEDGER) as f:
        return [json.loads(line) for line in f if line.strip()]


def spent(goal=None):
    """Money committed so far: real cost where known, the reservation otherwise."""
    open_reservations, total = {}, 0.0
    for e in _entries():
        if goal is not None and e["goal"] != goal:
            continue
        if e["kind"] == "reserve":
            open_reservations[e["id"]] = e["usd"]
        elif e["kind"] == "settle":
            open_reservations.pop(e["id"], None)
            total += e["usd"]
    return total + sum(open_reservations.values())


def spent_on_key(key):
    """Like spent(), for one key. Entries written before key 2 existed belong to key 1."""
    owner, open_reservations, total = {}, {}, 0.0
    for e in _entries():
        if e["kind"] == "reserve":
            owner[e["id"]] = e.get("key", "1")
            if owner[e["id"]] == key:
                open_reservations[e["id"]] = e["usd"]
        elif e["kind"] == "settle":
            open_reservations.pop(e["id"], None)
            if owner.get(e["id"], "1") == key:
                total += e["usd"]
    return total + sum(open_reservations.values())


def key_of(rid):
    for e in _entries():
        if e["kind"] == "reserve" and e["id"] == rid:
            return e.get("key", "1")
    return "1"


def remaining(goal):
    accounts = sum(max(0.0, cap - spent_on_key(k)) for k, cap in ACCOUNT_CAPS.items())
    return min(GOAL_CAPS[goal] - spent(goal), PROJECT_CAP - spent(), accounts)


def _append(entry):
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    with open(LEDGER, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write(json.dumps(entry) + "\n")
        f.flush()
        os.fsync(f.fileno())
        fcntl.flock(f, fcntl.LOCK_UN)


def reserve(goal, model, max_input_tokens, max_output_tokens, batch=False, note=""):
    """Claim the worst-case cost of a call. Raises OverBudget instead of overspending."""
    if goal not in GOAL_CAPS:
        raise KeyError(f"unknown goal {goal!r}")
    worst = cost(model, max_input_tokens, max_output_tokens, batch=batch)
    with _RESERVE_LOCK:
        left = min(GOAL_CAPS[goal] - spent(goal), PROJECT_CAP - spent())
        if worst > left:
            raise OverBudget(f"{goal}: call could cost ${worst:.2f} but only ${left:.2f} is left")
        key = next((k for k, cap in ACCOUNT_CAPS.items() if spent_on_key(k) + worst <= cap), None)
        if key is None:
            raise OverBudget(f"{goal}: no key has ${worst:.2f} left under its account cap")
        rid = f"{int(time.time() * 1000)}-{os.getpid()}-{os.urandom(3).hex()}"
        _append({"kind": "reserve", "id": rid, "goal": goal, "model": model, "usd": round(worst, 6),
                 "batch": batch, "note": note, "t": time.time(), "key": key})
    return rid


def settle(rid, goal, model, usage, batch=False):
    """Replace a reservation with the real cost, read from the API's usage object."""
    get = (lambda k: usage.get(k, 0) or 0) if isinstance(usage, dict) else (lambda k: getattr(usage, k, 0) or 0)
    usd = cost(model, get("input_tokens"), get("output_tokens"),
               get("cache_read_input_tokens"), get("cache_creation_input_tokens"), batch=batch)
    _append({"kind": "settle", "id": rid, "goal": goal, "model": model, "usd": round(usd, 6),
             "in": get("input_tokens"), "out": get("output_tokens"), "t": time.time()})
    return usd


def release(rid, goal):
    """A call that failed before the model ran costs nothing."""
    _append({"kind": "settle", "id": rid, "goal": goal, "model": "-", "usd": 0.0, "t": time.time()})


def report():
    lines = [f"{'goal':10s} {'spent':>9s} {'cap':>8s} {'left':>8s}"]
    for goal, cap in GOAL_CAPS.items():
        s = spent(goal)
        lines.append(f"{goal:10s} {s:9.2f} {cap:8.2f} {cap - s:8.2f}")
    lines.append(f"{'TOTAL':10s} {spent():9.2f} {PROJECT_CAP:8.2f} {PROJECT_CAP - spent():8.2f}")
    for k, cap in ACCOUNT_CAPS.items():
        lines.append(f"key {k}: ${cap - spent_on_key(k):.2f} left of ${cap:.2f}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(report())
