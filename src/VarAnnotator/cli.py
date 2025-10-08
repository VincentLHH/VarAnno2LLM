"""
VarAnnotator CLI: Variant annotator based on VEP
"""
import logging
import sys
from VarAnnotator.utils.logging import setup_logger
from VarAnnotator.preprocess.parser import create_parser
from .version import __version__

setup_logger()
logger = logging.getLogger(__name__)

def main():
    parser = create_parser()
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    setup_logger(level=log_level)
    logger.info(f"VarAnnotator v{__version__} start...")
    logger.debug(f"params: {args}")

    if args.command == "pre":
        from VarAnnotator.core import run_preprocessor
        success = run_preprocessor(args)
    elif args.command == "anno":
        from VarAnnotator.core import run_annotator
        success = run_annotator(args)
    elif args.command == "prompt":
        from VarAnnotator.core import run_prompter
        success = run_prompter(args)
    else:
        parser.print_help()
        exit(1)

    if success:
        sys.exit(0)
    else:
        sys.exit(1)


if __name__ == "__main__":
    main()