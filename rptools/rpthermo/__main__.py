#!/usr/bin/env python

from rptools.rpthermo  import runThermo, build_args_parser
import logging

def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    # Create logger
    logger = creage_logger('rptools - rpThermo')

    runThermo(,
            args.input_sbml,
            args.output_sink,
            args.remove_dead_end,
            args.compartment_id,
            logger=logger)


# used to initialise and download the data for equilibrator
def _init():
    from equilibrator_api import ComponentContribution
    cc = ComponentContribution()


def creage_logger(name):
    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(getattr(logging, args.log.upper()))
    return logger

if __name__ == '__main__':
    _init()
    _cli()
