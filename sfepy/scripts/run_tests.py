#!/usr/bin/env python
"""
Run SfePy tests. All arguments are passed to pytest.
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def main():
    parser = ArgumentParser(description=__doc__,
                            usage='sfepy-test [-h] [pytest options]',
                            formatter_class=RawDescriptionHelpFormatter)
    _, pytest_options = parser.parse_known_args()

    import sfepy
    sfepy.test(*pytest_options)

if __name__ == '__main__':
    main()
