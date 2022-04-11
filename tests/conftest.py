import pytest

@pytest.fixture(scope='session')
def output_dir(tmpdir_factory):
    """
    Output directory for tests.
    """
    return tmpdir_factory.mktemp('output')
