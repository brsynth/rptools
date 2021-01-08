#!/usr/bin/env python

from rptools.rpextractsink  import genSink, build_args_parser
from rptools.rplibs         import rpCache
import logging

def _cli():
    parser = build_args_parser()
    args  = parser.parse_args()

    # Create logger
    logger = logging.getLogger('rpExtractSink')
    logger.setLevel(getattr(logging, args.log.upper()))
    logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')

    rpcache = rpCache('file', ['cid_strc'], logger=logger)
    genSink(rpcache,
            args.input_sbml,
            args.output_sink,
            args.remove_dead_end,
            args.compartment_id,
            logger=logger)


if __name__ == '__main__':
    _cli()
