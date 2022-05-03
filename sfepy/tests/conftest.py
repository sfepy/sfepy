import pytest

def pytest_configure(config):
    config.addinivalue_line(
        'markers',
        'slow: Mark tests as slow (deselect with \'-m "not slow"\'.',
    )

def pytest_addoption(parser):
    parser.addoption('--output-dir', action='store', default=None)

@pytest.fixture(scope='session')
def output_dir(request, tmpdir_factory):
    """
    Output directory for tests.
    """
    output_dir = request.config.getoption('output_dir')
    if output_dir is not None:
        import os
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if not os.path.isdir(output_dir):
            raise IOError(f'cannot create directory "{output_dir}"!')

        return output_dir

    else:
        return tmpdir_factory.mktemp('output')
