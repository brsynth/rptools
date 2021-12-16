from os import (
    path as os_path,
    mkdir as os_mkdir
)
from logging import (
    StreamHandler,
    Logger,
    getLogger,
)
from colored import fg, bg, attr
from rr_cache import rrCache
from rptools import build_args_parser
from rptools.rpcompletion import rp_completion
from rptools.rpcompletion.Args import add_arguments


def _cli():
    parser = build_args_parser(
        prog='rpcompletion',
        description='Parse RP2 pathways to generate rpSBML collection of unique and complete (cofactors) pathways',
        m_add_args=add_arguments
    )
    args  = parser.parse_args()

    from rptools.__main__ import init
    logger = init(parser, args)


    check_args(
        args.max_subpaths_filter,
        args.outdir,
        logger
    )

    cache = rrCache(
        attrs=[
            'rr_reactions',
            'template_reactions',
            'cid_strc',
            'deprecatedCompID_compid',
        ]
    )

    pathways = rp_completion(
        rp2_metnet=args.rp2_metnet,
        sink=args.sink,
        rp2paths_compounds=args.rp2paths_compounds,
        rp2paths_pathways=args.rp2paths_pathways,
        cache=cache,
        upper_flux_bound=int(args.upper_flux_bound),
        lower_flux_bound=int(args.lower_flux_bound),
        max_subpaths_filter=args.max_subpaths_filter,
        logger=logger
    )

    # WRITE OUT
    if not os_path.exists(args.outdir):
        os_mkdir(args.outdir)
    # Write out selected pathways
    for pathway in pathways:
        pathway.to_rpSBML().write_to_file(
            os_path.join(
                args.outdir,
                'rp_'+pathway.get_id()
            ) + '.xml'
        )
    # for pathway_id, sub_pathways in pathways.items():
    #     for sub_pathway in sub_pathways[-args.max_subpaths_filter:]:
    #         sub_pathway.to_rpSBML().write_to_file(
    #             os_path.join(
    #                 args.outdir,
    #                 'rp_'+sub_pathway.get_id()
    #             ) + '.xml'
    #         )
    StreamHandler.terminator = ""
    logger.info(
        '{color}{typo}Results are stored in {rst}'.format(
            color=fg('white'),
            typo=attr('bold'),
            rst=attr('reset')
        )
    )
    StreamHandler.terminator = "\n"
    logger.info(
        '{color}{outdir}\n'.format(
            color=fg('grey_70'),
            outdir=args.outdir
        )
    )


def check_args(
    max_subpaths_filter: int,
    outdir: str,
    logger: Logger = getLogger(__name__)
):
    # out_format = out_format.upper()
    # if out_format not in FORMATS.keys():
    #     raise ValueError(
    #         'Output format {format} is not recognized (choices: {formats})'.format(
    #             format='\''+out_format+'\'',
    #             formats=', '.join(['\''+format+'\'' for format in FORMATS.keys()])
    #         )
    #     )
    #     exit(-1)

    if max_subpaths_filter < 0:
        raise ValueError('Max number of subpaths cannot be less than 0: '+str(max_subpaths_filter))

    if os_path.exists(outdir) and os_path.isfile(outdir):
        logger.error('Outdir name '+outdir+' already exists and is actually file. Stopping the process...')
        exit(-1)


if __name__ == '__main__':
    _cli()
