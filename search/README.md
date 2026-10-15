# The search system

How the results were found. Each target got a *chain*: a sequence of turns. In each turn:

1. A search program is written for the target, given the problem, the notes from earlier turns, the
   best program so far, and what the last program printed.
2. The program runs on a 64-core cloud machine for up to 40 minutes, without network access, as a
   separate user (`sandbox.py`).
3. Its answer is judged only by the target's own exact checker, never by the program itself.

Turns accumulate: a program can be an experiment whose printed measurements guide the next turn.
Programs were written with Anthropic's Claude Opus through the API; every call was charged to a
fixed budget before it was sent (`budget.py`, `llm.py`).

These files are reference copies from the full project: `chain.py` (the turn loop), `sandbox.py`
(local or remote isolated runs), `llm.py` and `budget.py` (API calls under a hard spending cap).
