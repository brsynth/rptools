#!/usr/bin/env python

from rptools.rpthermo import runThermo, build_args_parser
import logging

def _cli():
    parser = build_args_parser()
    args   = parser.parse_args()

    # Create logger
    logger = creage_logger('rptools - rpThermo', getattr(logging, args.log.upper()))

    runThermo(args.rpsbml_infile, args.rpsbml_outfile,
              args.pathway_id,
              args.ph, args.ionic_strength, args.pMg, args.temp_k,
              logger=logger)


def creage_logger(name=__name__, log_level='ERROR'):
    """
    Create a logger with name and log_level.

    Parameters
    ----------
    name : str
        A string containing the name that the logger will print out

    log_level : str
        A string containing the verbosity of the logger

    Returns
    -------
    Logger
        The logger object.

    """
    logger    = logging.getLogger(name)
    handler   = logging.StreamHandler()
    formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(log_level)

    return logger

if __name__ == '__main__':
    _cli()
