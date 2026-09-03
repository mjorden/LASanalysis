"""Thin wrapper kept for the Pages workflow; the builder lives in lasanalysis.site."""

import sys

from lasanalysis.site import main

if __name__ == "__main__":
    sys.exit(main())
