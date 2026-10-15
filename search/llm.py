"""The only place this project talks to the Claude API.

One model call at a time, each charged to a goal through common.budget. There
is deliberately no agent fan-out here: a call is a single request with a known
worst-case price, reserved before it is sent.

Needs `pip install anthropic` and ANTHROPIC_API_KEY in the environment.
"""
from . import budget

DEFAULT_MODEL = "claude-opus-5-5"


def _load_env():
    """Read ANTHROPIC_API_KEY from an untracked .env file at the repo root, if it is not already set."""
    import os
    path = os.path.join(budget.ROOT, ".env")
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                key, sep, value = line.strip().partition("=")
                if sep and key in ("ANTHROPIC_API_KEY", "ANTHROPIC_API_KEY_2") and key not in os.environ:
                    os.environ[key] = value.strip().strip('"').strip("'")


def _client(key="1"):
    """A client for key 1 (ANTHROPIC_API_KEY) or key 2 (ANTHROPIC_API_KEY_2)."""
    import anthropic  # imported lazily so the offline parts of the repo run without the SDK
    import os
    _load_env()
    return anthropic.Anthropic(api_key=os.environ["ANTHROPIC_API_KEY" if key == "1" else "ANTHROPIC_API_KEY_2"])


def ask(goal, prompt, *, system=None, model=DEFAULT_MODEL, effort="high", max_tokens=32000, note=""):
    """Send one prompt and return (text, usd). Raises budget.OverBudget before overspending.

    Thinking is adaptive on every model used here and is billed as output, so
    `max_tokens` is the worst-case output bill for the call.
    """
    import anthropic

    est_input = (len(prompt) + len(system or "")) // 3 + 200
    rid = budget.reserve(goal, model, est_input, max_tokens, note=note)
    kwargs = dict(
        model=model,
        max_tokens=max_tokens,
        thinking={"type": "adaptive"},
        output_config={"effort": effort},
        messages=[{"role": "user", "content": prompt}],
    )
    if system:
        # The system prompt is the stable prefix across a search run: cache it.
        kwargs["system"] = [{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}]
    import time

    message = None
    for attempt in range(6):
        try:
            with _client(budget.key_of(rid)).messages.stream(**kwargs) as stream:
                message = stream.get_final_message()
            break
        except Exception as e:
            # A dropped stream (connection reset mid-reply) surfaces as the HTTP library's own ReadError.
            if not isinstance(e, (anthropic.APIConnectionError, anthropic.InternalServerError, anthropic.RateLimitError)) \
                    and type(e).__name__ not in ("ReadError", "RemoteProtocolError"):
                if isinstance(e, anthropic.APIStatusError):
                    budget.release(rid, goal)
                raise
            # A network blip or a busy server: wait and try again before giving up on a long run.
            if attempt == 5:
                budget.release(rid, goal)
                raise
            time.sleep(30 * (attempt + 1))
    usd = budget.settle(rid, goal, model, message.usage)
    if message.stop_reason == "refusal":
        return "", usd
    text = "".join(block.text for block in message.content if block.type == "text")
    return text, usd


def submit_batch(goal, prompts, *, system=None, model=DEFAULT_MODEL, effort="medium", max_tokens=16000, note="",
                 shared_prefix=None):
    """Queue many independent prompts at half price. `prompts` maps custom_id -> prompt.

    Returns (batch_id, reservation_id). Collect with `collect_batch` once it has ended.
    """
    from anthropic.types.message_create_params import MessageCreateParamsNonStreaming
    from anthropic.types.messages.batch_create_params import Request

    est_input = (sum(len(p) for p in prompts.values()) // 3
                 + len(prompts) * (200 + (len(system or "") + len(shared_prefix or "")) // 3))
    rid = budget.reserve(goal, model, est_input, max_tokens * len(prompts), batch=True, note=note)
    requests = []
    for custom_id, prompt in prompts.items():
        params = dict(
            model=model,
            max_tokens=max_tokens,
            thinking={"type": "adaptive"},
            output_config={"effort": effort},
            messages=[{"role": "user", "content": prompt}],
        )
        if shared_prefix:
            # Every request in a generation starts with the same task and best programs. Marking that
            # part cacheable lets later requests read it at a tenth of the input price. Inside a batch
            # the hits are best-effort, so this trims the input bill; it does not remove it.
            params["messages"] = [{"role": "user", "content": [
                {"type": "text", "text": shared_prefix, "cache_control": {"type": "ephemeral"}},
                {"type": "text", "text": prompt},
            ]}]
        if system:
            params["system"] = [{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}]
        requests.append(Request(custom_id=custom_id, params=MessageCreateParamsNonStreaming(**params)))
    import time
    import anthropic

    for attempt in range(6):
        try:
            batch = _client(budget.key_of(rid)).messages.batches.create(requests=requests)
            return batch.id, rid
        except anthropic.APITimeoutError:
            # The request may have reached the server before the reply was lost. Repeating it could queue
            # a second, paid copy, so stop here and leave the hold on the ledger as the cautious figure.
            raise
        except anthropic.APIConnectionError:
            # Never connected, so nothing was queued: safe to try again.
            if attempt == 5:
                budget.release(rid, goal)
                raise
            time.sleep(30 * (attempt + 1))
        except anthropic.APIStatusError:
            budget.release(rid, goal)
            raise


def collect_batch(goal, batch_id, rid, *, model=DEFAULT_MODEL, settle=True):
    """Return {custom_id: text} for an ended batch, or None if it is still running."""
    client = _client(budget.key_of(rid) if rid else "1")
    if client.messages.batches.retrieve(batch_id).processing_status != "ended":
        return None
    texts = {}
    totals = {"input_tokens": 0, "output_tokens": 0, "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0}
    for result in client.messages.batches.results(batch_id):   # arrives in any order: key by custom_id
        if result.result.type != "succeeded":
            texts[result.custom_id] = None
            continue
        message = result.result.message
        for key in totals:
            totals[key] += getattr(message.usage, key, 0) or 0
        texts[result.custom_id] = "".join(b.text for b in message.content if b.type == "text")
    if settle:               # re-reading a finished batch is free and must not be charged twice
        budget.settle(rid, goal, model, totals, batch=True)
    return texts
