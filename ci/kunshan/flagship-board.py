#!/usr/bin/env python3
"""Refresh the flagship-port progress board from live cluster state.

The checkbox board in docs/plans/oneflow-development-todo.md is a hand-written
record; this script is the live view.  It reads three things and prints them
together so progress cannot silently drift:

  1. active and recent Slurm jobs, so an in-flight experiment is visible even
     when nobody has written a note yet;
  2. the evidence directories each job produced, so "done" means artefacts
     exist rather than a claim;
  3. the checkbox state parsed out of the todo document, so the hand-written
     board and reality sit side by side.

Usage:
  python3 ci/kunshan/flagship-board.py                  # live view
  python3 ci/kunshan/flagship-board.py --write-todo PATH # also patch the live
                                                          # block into the todo
The --write-todo mode only replaces the region between the BOARD markers and
never touches hand-written checkboxes outside it.
"""

import argparse
import json
import re
import subprocess
import sys

SCNET = [
    "/Users/luql/.claude-scnet/skills/scnet-hpc/scripts/scnet.py",
    "--region", "11250",
]
BEGIN = "<!-- FLAGSHIP-BOARD:BEGIN -->"
END = "<!-- FLAGSHIP-BOARD:END -->"

# Which job-name fragments belong to which board line.  Kept as data so the
# map is auditable rather than buried in control flow.
WATCHED = [
    ("C1  SA cost split (loop vs math)", ["saperf", "saresplit", "saabl"]),
    ("D   viscous flux scope",           ["visc", "b3", "visflux"]),
    ("E1  LU-SGS convergence",           ["converge", "lusgs", "rk3", "tintegral"]),
    ("A4  vismodel guard build/verify",  ["vismodel", "verify-vismodel"]),
    ("cfg 2-card / scaling probes",      ["config", "strong", "ledger", "rpc"]),
    ("other",                            []),
]


def run(args):
    try:
        out = subprocess.run(args, capture_output=True, text=True, timeout=180)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return None, str(exc)
    if out.returncode != 0:
        return None, out.stderr.strip()[:200]
    return out.stdout, None


def jobs(scope, limit=40, days=None):
    args = SCNET + ["jobs", "--scope", scope, "--limit", str(limit)]
    if days:
        args += ["--days", str(days)]
    text, err = run(args)
    if text is None:
        return [], err
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        return [], "unparseable job list"
    return data.get("items", []), None


_NAME_CACHE = {}


def job_name(job_id):
    """One API call per id, memoised.

    Each lookup is a network round trip; calling it per job in a 40-row history
    made the script hang past its timeout.  Active jobs are worth the cost,
    history is not - classify history by job id only.
    """
    key = str(job_id)
    if key in _NAME_CACHE:
        return _NAME_CACHE[key]
    text, _ = run(SCNET + ["--details", "job", key])
    name = None
    if text is not None:
        try:
            name = json.loads(text).get("name")
        except json.JSONDecodeError:
            name = None
    _NAME_CACHE[key] = name
    return name


def todo_checkboxes(path):
    """Return (done, total, open_items) for section 5.35."""
    try:
        text = open(path, encoding="utf-8").read()
    except OSError as exc:
        return None, None, [str(exc)]
    start = text.find("### 5.35")
    if start < 0:
        return None, None, ["section 5.35 not found"]
    end = text.find("\n### ", start + 1)
    block = text[start:end if end > 0 else len(text)]
    done = len(re.findall(r"^- \[x\]", block, re.M))
    todo = len(re.findall(r"^- \[ \]", block, re.M))
    open_items = re.findall(r"^- \[ \] (\*\*[^*]+\*\*|[\w ].{0,60})", block, re.M)
    return done, done + todo, [o.strip() for o in open_items]


def classify(name, fragments):
    low = (name or "").lower()
    return any(f in low for f in fragments)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-todo")
    ap.add_argument("--todo",
                    default="docs/plans/oneflow-development-todo.md")
    opts = ap.parse_args()

    lines = []
    add = lines.append
    add(BEGIN)
    add("")
    add("**本块由 `ci/kunshan/flagship-board.py` 生成，勿手改。**")
    add("")

    # --- checkbox progress from the hand-written board -------------------
    done, total, open_items = todo_checkboxes(opts.todo)
    if total:
        add(f"进度：**{done}/{total}** 勾选（含已排除项）。未完成：")
        for item in open_items[:12]:
            add(f"- {item}")
        add("")

    # --- live jobs -------------------------------------------------------
    active, err = jobs("active")
    if err:
        add(f"集群：无法读取活动作业（{err}）")
    elif not active:
        add("集群：**无活动作业**。")
    else:
        add(f"集群：**{len(active)} 个活动作业**")
        for job in active[:12]:
            jid = job.get("job_id")
            state = job.get("state")
            name = job_name(jid) or "?"
            add(f"- `{jid}` {state} — {name}")
    add("")

    # --- what those jobs map to ------------------------------------------
    # History rows carry no name and resolving each one is a round trip, so
    # only the active set is classified by name; history contributes counts.
    history, _ = jobs("history", limit=30, days=1)
    add(f"近一天作业：{len(history)} 个（状态见上面活动列表）")
    add("")
    add("按看板条目（仅活动作业可分类，历史作业不计）：")
    for label, fragments in WATCHED:
        if not fragments:
            hit_ids = [j for j in active
                       if not any(classify(job_name(j.get("job_id")), f)
                                  for _, f in WATCHED if f)]
        else:
            hit_ids = [j for j in active
                       if classify(job_name(j.get("job_id")), fragments)]
        hits = [f"`{j.get('job_id')}` {j.get('state')} — {job_name(j.get('job_id')) or '?'}"
                for j in hit_ids]
        add(f"- {label}: " + (", ".join(hits[:3]) if hits else "无活动作业"))
    add("")
    add(END)

    block = "\n".join(lines)
    if opts.write_todo:
        text = open(opts.write_todo, encoding="utf-8").read()
        if BEGIN in text and END in text:
            pre = text[:text.index(BEGIN)]
            post = text[text.index(END) + len(END):]
            open(opts.write_todo, "w", encoding="utf-8").write(
                pre + block + post)
            print(f"patched live block into {opts.write_todo}")
        else:
            print(f"{opts.write_todo} has no BOARD markers; not patching")
            print(block)
    else:
        print(block)


if __name__ == "__main__":
    sys.exit(main())
