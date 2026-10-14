"""Run a model-written Python program without trusting it.

The program gets a time limit, no network, and can write only inside its own
scratch directory (macOS `sandbox-exec`). It reports its answer by printing one
line of JSON as the last line of stdout; everything else it prints is ignored.
The answer is then checked by our own verifier, never taken on trust.
"""
import json
import os
import subprocess
import sys
import tempfile

PROFILE = """(version 1)
(allow default)
(deny network*)
(deny file-write*)
(allow file-write* (subpath "{scratch}") (literal "/dev/null") (literal "/dev/stdout") (literal "/dev/stderr"))
"""


def run(source, seconds, memory_note=""):
    """Returns (answer, log). `answer` is the parsed last JSON line, or None."""
    with tempfile.TemporaryDirectory(prefix="moonshot-") as scratch:
        scratch = os.path.realpath(scratch)
        path = os.path.join(scratch, "candidate.py")
        with open(path, "w") as f:
            f.write(source)
        cmd = ["sandbox-exec", "-p", PROFILE.format(scratch=scratch), sys.executable, "-I", path]
        try:
            proc = subprocess.run(cmd, cwd=scratch, capture_output=True, text=True, timeout=seconds)
        except subprocess.TimeoutExpired as e:
            out = e.stdout.decode() if isinstance(e.stdout, bytes) else (e.stdout or "")
            return _last_json(out), f"timed out after {seconds}s"
        log = (proc.stderr or "")[-2000:]
        if proc.returncode != 0:
            return _last_json(proc.stdout), f"exit {proc.returncode}: {log}"
        return _last_json(proc.stdout), log


def _last_json(stdout):
    for line in reversed((stdout or "").strip().splitlines()):
        line = line.strip()
        if line.startswith(("[", "{")):
            try:
                return json.loads(line)
            except json.JSONDecodeError:
                continue
    return None


def run_full(source, seconds):
    """Like `run`, but also returns what the program printed: (answer, stdout tail, stderr tail).

    For searches that work in turns and read their own experiments back. C may be compiled inside
    the scratch directory, so the compiler's temporary files are pointed there too.
    """
    if os.environ.get("MOONSHOTS_REMOTE"):
        return _run_remote(source, seconds)
    with tempfile.TemporaryDirectory(prefix="moonshot-") as scratch:
        scratch = os.path.realpath(scratch)
        path = os.path.join(scratch, "candidate.py")
        with open(path, "w") as f:
            f.write(source)
        cmd = ["sandbox-exec", "-p", PROFILE.format(scratch=scratch), sys.executable, "-I", path]
        env = dict(os.environ, TMPDIR=scratch)
        env.pop("ANTHROPIC_API_KEY", None)
        try:
            proc = subprocess.run(cmd, cwd=scratch, capture_output=True, text=True, timeout=seconds, env=env)
            out, err = proc.stdout or "", proc.stderr or ""
            if proc.returncode != 0:
                err += f"\n[exit code {proc.returncode}]"
        except subprocess.TimeoutExpired as e:
            text = lambda b: b.decode(errors="replace") if isinstance(b, bytes) else (b or "")
            out, err = text(e.stdout), text(e.stderr) + f"\n[killed at the {seconds}s limit]"
        return _last_json(out), out[-3000:], err[-1500:]


def _run_remote(source, seconds):
    """Run the program on the rented machine named by MOONSHOTS_REMOTE (user@host) instead of here.

    There it runs as the user `runner`, whose outgoing network traffic the firewall rejects, in a
    directory of its own that is deleted afterwards. The API key never leaves this computer.
    """
    host = os.environ["MOONSHOTS_REMOTE"]
    ssh = ["ssh", "-i", os.environ["MOONSHOTS_REMOTE_KEY"], "-o", "ServerAliveInterval=30", "-o", "ConnectTimeout=30",
           "-o", "BatchMode=yes", host]
    job = "/work/" + os.urandom(6).hex()
    try:
        subprocess.run(ssh + [f"mkdir -p {job} && chmod 777 {job} && cat > {job}/candidate.py"], input=source, text=True,
                       timeout=120, check=True, capture_output=True)
        proc = subprocess.run(ssh + [f"cd {job} && sudo -u runner env TMPDIR={job} HOME={job} timeout -k 5 {int(seconds)} "
                                     f"python3 -I candidate.py; code=$?; sudo pkill -9 -u runner -f {job} ; exit $code"],
                              capture_output=True, text=True, timeout=seconds + 180)
        out, err = proc.stdout or "", proc.stderr or ""
        if proc.returncode == 124:
            err += f"\n[stopped at the {seconds}s limit]"
        elif proc.returncode != 0:
            err += f"\n[exit code {proc.returncode}]"
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as e:
        out, err = "", f"[could not run on the remote machine: {type(e).__name__}]"
    subprocess.run(ssh + [f"sudo rm -rf {job}"], capture_output=True, timeout=120)
    return _last_json(out), out[-3000:], err[-1500:]
