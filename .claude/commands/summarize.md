---
name: summarize
description: provides a summary of accomplishments and struggles occurred in current session
---

# Session Summary

## Purpose
Provide a concise 3-sentence summary of what has been done in the current session to inform a supervisor of the user's accomplishment and struggles.

## Process
1. Look at what has been done in the current session — files edited, commands run, questions asked.
2. Check the actual work artifacts, not just the chat:
   - If a notebook is open or was worked on this session, read its code cells and markdown-answer cells directly. Compare each against its placeholder text (e.g. "Your final ... script here", "(Write your answers here)") to tell filled-in work apart from untouched cells.
   - Run `git status` and `git diff` to see which files were actually modified or created — this can surface changes made outside the chat.
   - Scan cell outputs for tracebacks, error messages, or failed assertions — these count as struggles even if the intern never mentioned them in chat.
3. Write the summary as three separate bullet points, one sentence each, in this format:
   - Accomplished: <what was accomplished — cite specific exercises/files completed, based on the artifact check above, not just what was discussed>
   - Broke/surprising: <anything that broke or was surprising — including errors found in cell outputs, not only ones raised in conversation>
   - Still open: <what is still open or unfinished — explicitly name any exercise/cell still containing placeholder text>
4. Print the three bullet points to the terminal so the intern can copy them into a note or share them with their supervisor
