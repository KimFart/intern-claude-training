---
name: log
description: End-of-session notebook. Run at the end of every work session to record what you did, what broke, and what you learned. Builds a running lab notebook in log/session-log.md.
---

# Session Log

## Purpose

Record what happened this session. This is your lab notebook for code. Future you — and future interns — will thank you.

## Process

1. Read the existing `log/session-log.md` if it exists to understand the running context.

2. Ask the user three questions (one at a time if they haven't answered already):
   - What did you accomplish this session?
   - What broke or didn't work as expected?
   - What did you learn or understand better now?

3. Append a new entry to `log/session-log.md` in this format:

```markdown
## Session — YYYY-MM-DD

### Done
- [bullet points of what was accomplished]

### Broke / Struggled
- [bullet points of errors, failures, confusing things]

### Learned
- [bullet points of new understanding]

---
```

4. Create `log/` directory if it doesn't exist.

5. Confirm: "Logged. Your session notebook is at `log/session-log.md`."

## Notes

- Be concise. Each bullet should be one sentence.
- "Broke" entries are not failures — they are the most useful part of the log.
- Do not summarize or editorialize. Record what the user tells you.
