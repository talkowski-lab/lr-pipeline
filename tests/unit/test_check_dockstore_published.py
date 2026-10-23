import subprocess
import sys
import urllib.error

import pytest
from conftest import load_script

pub = load_script(".github/scripts/check_dockstore_published.py")


def http_error(code):
    return urllib.error.HTTPError("https://dockstore.org", code, "msg", None, None)


class Response:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


@pytest.fixture
def no_sleep(monkeypatch):
    monkeypatch.setattr(pub.time, "sleep", lambda _: None)


def fake_urlopen(monkeypatch, outcomes):
    """Replace urlopen with a queue of responses or exceptions, recording the requested URLs."""
    urls = []
    queue = list(outcomes)

    def urlopen(url, timeout):
        urls.append(url)
        outcome = queue.pop(0)
        if isinstance(outcome, Exception):
            raise outcome
        return outcome

    monkeypatch.setattr(pub.urllib.request, "urlopen", urlopen)
    return urls


def test_is_published_true_on_200_and_url_is_encoded(monkeypatch):
    urls = fake_urlopen(monkeypatch, [Response()])
    assert pub.is_published("TRGT") is True
    assert urls == ["https://dockstore.org/api/workflows/path/workflow/github.com%2Ftalkowski-lab%2Flr-pipeline%2FTRGT/published"]


def test_is_published_false_on_404(monkeypatch):
    fake_urlopen(monkeypatch, [http_error(404)])
    assert pub.is_published("TRGT") is False


def test_is_published_retries_server_errors_and_timeouts(monkeypatch, no_sleep):
    urls = fake_urlopen(monkeypatch, [http_error(503), TimeoutError(), Response()])
    assert pub.is_published("TRGT") is True
    assert len(urls) == 3


def test_is_published_gives_up_after_repeated_failures(monkeypatch, no_sleep):
    fake_urlopen(monkeypatch, [http_error(500)] * pub.ATTEMPTS)
    with pytest.raises(RuntimeError):
        pub.is_published("TRGT")


def test_find_added_names_parses_diff(monkeypatch):
    diff = "+  - subclass: WDL\n+    name: NewOne\n-    name: Gone\n     name: Unchanged\n+    name: AnotherNew\n"
    result = subprocess.CompletedProcess([], 0, stdout=diff)
    monkeypatch.setattr(pub.subprocess, "run", lambda *a, **k: result)
    assert pub.find_added_names("HEAD~1") == {"NewOne", "AnotherNew"}


@pytest.fixture
def run(monkeypatch, capsys):
    monkeypatch.setattr(pub, "load_workflow_names", lambda: ["A", "B", "C"])

    def go(published, *args, added=()):
        monkeypatch.setattr(pub, "is_published", lambda name: name in published)
        monkeypatch.setattr(pub, "find_added_names", lambda ref: set(added))
        monkeypatch.setattr(sys, "argv", ["check_dockstore_published.py", *args])
        try:
            pub.main()
            code = 0
        except SystemExit as exc:
            code = exc.code
        return code, capsys.readouterr().out

    return go


def test_main_passes_when_all_published(run):
    code, out = run({"A", "B", "C"})
    assert code == 0
    assert "passed" in out


def test_main_fails_listing_unpublished(run):
    code, out = run({"A"})
    assert code == 1
    assert "'B'" in out and "'C'" in out and "'A'" not in out


def test_main_added_since_downgrades_new_entries_to_warnings(run):
    code, out = run({"A"}, "--added-since", "abc", added={"B", "C"})
    assert code == 0
    assert "warning: 'B'" in out and "warning: 'C'" in out


def test_main_added_since_still_fails_on_older_unpublished(run):
    code, out = run({"A"}, "--added-since", "abc", added={"B"})
    assert code == 1
    assert "warning: 'B'" in out
    assert "- 'C' is not published" in out
