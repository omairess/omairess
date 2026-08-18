# Downloads

## BSRT-webapp.zip — the browser version, ready to run

**To download from GitHub:** click `BSRT-webapp.zip` above, then the **Download**
button (GitHub shows binary files as "view raw" rather than rendering them).

Unzip it anywhere and double-click `BSRT/index.html`. It runs offline in Chrome,
Edge or Firefox — no installation, no server, no account. `OPEN_ME_FIRST.txt`
inside covers the first five minutes; `README.md` is the full documentation.

This is a **copy of `/BSRT` at the commit that last touched it**, bundled for
people who should not have to use git. If you are comfortable with git, cloning
the repository and opening `BSRT/index.html` gives you exactly the same app and
stays current automatically.

### Keeping it current

The zip does not rebuild itself. When `/BSRT` changes, this file is stale until
someone regenerates it:

```bash
cd /path/to/repo
rm -rf /tmp/z && mkdir -p /tmp/z && cp -r BSRT /tmp/z/
cp downloads/OPEN_ME_FIRST.txt /tmp/z/BSRT/ 2>/dev/null || true
(cd /tmp/z && zip -rq BSRT-webapp.zip BSRT)
mv /tmp/z/BSRT-webapp.zip downloads/
```

## BSRT-psychopy.zip — the PsychoPy version

For labs that already run PsychoPy, or are willing to install it. It is the
most accurate of the three builds: stimulus onsets come from the vsync
timestamp of the frame that carried them, and reaction times from the
keyboard's own hardware stamp.

**Download it the same way**, unzip it, and open `START_HERE.txt` — a guide
written for someone who has never used PsychoPy. It walks through installing
the standalone application, running `check_setup.py` to confirm the timing on
that machine is trustworthy, and then running the task.

**Run `check_setup.py` before any participant.** It takes a minute and tells
you whether the screen clock and the keyboard clock agree, which is the one
thing that would otherwise be wrong invisibly.

This build has not yet been run against a real PsychoPy installation — see the
note at the end of `START_HERE.txt`.

## The desktop version is not here

macOS and Windows builds cannot be produced on this machine — packaging has to
run on the target operating system. Get them from the **Actions** tab
(*Build BSRT desktop apps* → *Run workflow*), or from the Releases page if a
version tag has been pushed. See `BSRT-desktop/README.md`.
