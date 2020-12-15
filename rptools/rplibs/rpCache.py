from os         import path as os_path
from os         import mkdir as os_mkdir
from rdkit.Chem import MolFromSmiles, MolFromInchi, MolToSmiles, MolToInchi, MolToInchiKey
from csv        import DictReader as csv_DictReader
from csv        import reader as csv_reader
from json       import dump as json_dump
from json       import load as json_load
from gzip       import open as gzip_open
from re         import findall as re_findall
from time       import time as time_time
from brs_utils  import print_OK, print_FAILED, download
from requests   import exceptions as r_exceptions
from redis      import StrictRedis
from credisdict import CRedisDict, wait_for_redis
from hashlib    import sha512
from pathlib    import Path
from colored    import attr as c_attr
import logging


#######################################################
################### rpCache  ##########################
#######################################################


class rpCache:
    """Class to generate the cache

    Contains all the functions that parse different files, used to calculate the thermodynamics and the FBA of the the other steps. These should be called only when the files have changes

    """
    # _input_cache_url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'
    _cache_url       = 'https://gitlab.com/breakthewall/rpcache-data/-/raw/master/'

    # static attribues
    _convertMNXM = {
        'MNXM162231': 'MNXM6',
        'MNXM84':     'MNXM15',
        'MNXM96410':  'MNXM14',
        'MNXM114062': 'MNXM3',
        'MNXM145523': 'MNXM57',
        'MNXM57425':  'MNXM9',
        'MNXM137':    'MNXM588022'
        }

    # name: sha512sum
    _input_cache_files = {
            'chem_xref.tsv.gz':    'e558110990dcc75af943863790dc55360fd2d40ecb17d02335377671e80f0ab3738fd556acb340e03e48dd1afdec3eece1e92df1e18bc24e7445f24f778a10da',
            'reac_xref.tsv.gz':    '48b991cf4a9c2ca573d395cf35c378881ed79e87772827647bfab2f6345499698664e07195ec10b342fc0164304dbd2363cccff1a1182225e6afebce3c16448b',
            'compounds.tsv.gz':    '719716bb880257bd014e045c03eb8dd12e2bbeba3aa52e38e9632ce605817b9dc09530e81fadd25542c0a439bdb81e1dfbd3a38f35b30b061845d1a880dbfe01',
            'chem_prop.tsv.gz':    'f2d220d1f0425e5e47f01e7deccfa46b60094d43b9f62b191ffb0fab8c00ef79e87c3b71d10bdcd26020608094f24884f51b3ebc3d7d3c9a6d594c6eaa324c66',
            'retrorules_rr02_flat_all.tsv.gz':   '890bdd24042c0192b5538964d775feefcb6cff9ad5f35690bfbfc5ae09334dd19df6828cdfc7f57a2018e090571517122b99d8760128052af898c638ae667e24',
            'comp_xref.tsv.gz':    '913a827f3645fda1699676ae6c32b9d7a8debae97ce7b0c386d8447f4eee5aa721d31bfb856d4092b3d5e987a8f19a6fe4bd28ddf1c5df5f85e71c3625bd1d81',
            'rxn_recipes.tsv.gz':  'dc0624f5ed7ab0b691d9a6ba02571a5cf334cfdb3109e78c98708e31574c46aeac2a97e9433788d80490ff80337679ccfd706cbb8e71a11cdc6122573bb69b0f'
            }

    # Attributes with dependencies (other attributes + input_cache files)
    __attributes = {
            'deprecatedCID_cid': {'attr_deps': [],
                                  'file_deps': ['chem_xref.tsv.gz']},

            'deprecatedRID_rid': {'attr_deps': [],
                                  'file_deps': []},

            'cid_strc': {'attr_deps': ['deprecatedCID_cid'],
                         'file_deps': ['compounds.tsv.gz', 'chem_prop.tsv.gz']},

            'cid_name': {'attr_deps': ['deprecatedCID_cid'],
                         'file_deps': ['compounds.tsv.gz', 'chem_prop.tsv.gz']},

            'cid_xref': {'attr_deps': ['deprecatedCID_cid'],
                         'file_deps': []},

            'chebi_cid': {'attr_deps': ['cid_xref'],
                          'file_deps': []},

            'rr_reactions': {'attr_deps': ['deprecatedCID_cid', 'deprecatedRID_rid'],
                             'file_deps': ['retrorules_rr02_flat_all.tsv.gz']},

            'inchikey_cid': {'attr_deps': ['cid_strc'],
                             'file_deps': []},

            'comp_xref': {'attr_deps': [],
                          'file_deps': ['comp_xref.tsv.gz']},

            'deprecatedCompID_compid': {'attr_deps': [],
                                        'file_deps': ['comp_xref.tsv.gz']},

            'rr_full_reactions': {'attr_deps': ['deprecatedCID_cid', 'deprecatedRID_rid'],
                                  'file_deps': ['rxn_recipes.tsv.gz']},
    }

    _attributes = list(__attributes.keys())


    # name: sha512sum
    _cache_files = {
            _attributes[0]+'.json.gz': '698a3e83cf4f9206ea2644c9c35a9af53957838baaae6efb245d02b6b8d0ea8b25c75008e562b99ba3e0189e50ee47655376f2d0635f6206e0015f91f0e4bad8',
            _attributes[1]+'.json.gz': '51554c6f6ae99c6755da7496208b3feec30547bc4cf3007d9fd30f46fa4c0cc73bad5aeb743dca07e32711c4346504296bee776d135fb18e96c891a0086fc87e',
            _attributes[2]+'.json.gz': '0021ef63165d75ee6b8c209ccf14b8a1b8b7b263b4077f544729c47b5525f66511c3fa578fd2089201abb61693085b9912639e62f7b7481d06ad1f38bfc2dd8e',
            _attributes[3]+'.json.gz': '7d559cc7389c0cb2bd10f92e6e845bb5724be64d1624adc4e447111fc63599bb69396cd0cc3066a6bb19910c00e266c97e21b1254d9a6dc9da3a8b033603fcff',
            _attributes[4]+'.json.gz': '587d6c5206ee94e63af6d9eaf49fd5e2ca417308b3ece8a7f47e916c42376e2c8635a031ce26dc815cd7330f2323054a44d23951e416a9a29c5a9a2ab51e8953',
            _attributes[5]+'.json.gz': '8783aaa65a281c4a7ab3a82a6dc99620418ed2be4a739f46db8ee304fcb3536a78fed5a955e1c373a20c3e7d3673793157c792b4429ecb5c68ddaddb1a0f7de7',
            _attributes[6]+'.json.gz': '8007480fc607caf41f0f9a93beb66c7caa66c37a3d01a809f6b94bc0df469cec72091e8cc0fbabb3bd8775e9776b928ecda2779fc545c7e4b9e71c504f9510ce',
            _attributes[7]+'.json.gz': 'afc2ad3d31366a8f7fe1604fa49c190ade6d46bc8915f30bd20fdfdfc663c979bb10ca55ad10cadec6002a17add46639c70e7adf89cb66c57ed004fd3e4f0051',
            _attributes[8]+'.json.gz': '81c673fe1940e25a6a9722fd74b16bc30e1590db0c40810f541ad4ffba7ae04c01268b929d4bf944e84095a0c2a1d0079d1861bc1df3e8308fbb6b35e0aaf107',
            _attributes[9]+'.json.gz': '599e4de4935d2ba649c0b526d8aeef6f0e3bf0ed9ee20adad65cb86b078ac139e4cc9758945c2bb6da1c6840867239c5415cb5bceeb80164798ff627aac0a985',
            _attributes[10]+'.json.gz': '599e4de4935d2ba649c0b526d8aeef6f0e3bf0ed9ee20adad65cb86b078ac139e4cc9758945c2bb6da1c6840867239c5415cb5bceeb80164798ff627aac0a985'
            }



    _ext = '.json.gz'

    _pubchem_species = {}

    ## Cache constructor
    #
    # @param self The object pointer
    # @param db Mode of storing objects ('file' or 'redis')
    def __init__(self, db='file', attrs='', logger=None):

        if logger is None:
            # Create logger
            self.logger = logging.getLogger(__name__)
            self.logger.setLevel(getattr(logging, 'ERROR'))
            self.logger.formatter = logging.Formatter('%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s')
        else:
            self.logger = logger

        self.logger.info('Started instance of rpCache')

        self.store_mode = db
        rpCache._db_timeout = 10

        if attrs:
            self._attributes = attrs

        self.dirname = os_path.dirname(os_path.abspath( __file__ ))#+"/.."
        # input_cache
        self._input_cache_dir = self.dirname+'/input_cache/'
        # cache
        self._cache_dir = self.dirname+'/cache/'

        if self.store_mode!='file':
            self.redis = StrictRedis(host=self.store_mode, port=6379, db=0, decode_responses=True)
            if not wait_for_redis(self.redis, self._db_timeout):
                self.logger.critical("Database "+self.store_mode+" is not reachable")
                exit()
            for attr in self._attributes:
                setattr(self, attr, CRedisDict(attr, self.redis))
        else:
            for attr in self._attributes:
                setattr(self, attr, None)

        try:
            self._check_or_load_cache()
        except FileNotFoundError:
            print_FAILED()
            try:
                rpCache._check_or_download_cache_to_disk(self._cache_dir, self._attributes)
                self._check_or_load_cache()
            except (r_exceptions.RequestException,
                    r_exceptions.InvalidSchema,
                    r_exceptions.ConnectionError):
                print_FAILED()
                rpCache.generate_cache(self._cache_dir)
                self._check_or_load_cache()


    def get(self, attr):
        return getattr(self, attr)

    #####################################################
    ################# ERROR functions ###################
    #####################################################

    class Error(Exception):
        """Error function for the convertion of structures
        """
        pass


    class DepictionError(Error):
        """Error function for the convertion of structures
        """
        def __init__(self, message):
            """Constructor for the class

            :param message: The error handling message string

            :type message: str

            :rtype: None
            :return: None
            """
            #self.expression = expression
            self.message = message

    #url = 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/'
    #url = 'ftp://ftp.vital-it.ch/databases/metanetx/MNXref/3.2/'

    @staticmethod
    def generate_cache(outdir, logger=None):

        logger = logger or logging.getLogger(__name__)

        if outdir == '':
            outdir = 'cache'

        if not os_path.isdir(outdir):
            os_mkdir(outdir)
        outdir += '/'

        url = rpCache._cache_url

        # FETCH INPUT_CACHE FILES
        input_dir = 'input-'+os_path.basename(os_path.normpath(outdir))+'/'
        for file in rpCache._input_cache_files.keys():
            rpCache._download_input_cache(url, file, input_dir)

        # GENERATE CACHE FILES AND STORE THEM TO DISK
        deprecatedCID_cid  = rpCache._gen_deprecatedCID_cid(input_dir, outdir)
        cid_strc, cid_name = rpCache._gen_cid_strc_cid_name(input_dir, outdir, deprecatedCID_cid)
        rpCache._gen_inchikey_cid(input_dir, outdir, cid_strc)
        del cid_strc, cid_name
        cid_xref           = rpCache._gen_cid_xref(input_dir, outdir, deprecatedCID_cid)
        rpCache._gen_chebi_cid(input_dir, outdir, cid_xref)
        del cid_xref
        deprecatedRID_rid  = rpCache._gen_deprecatedRID_rid(input_dir, outdir)
        rpCache._gen_rr_reactions(input_dir, outdir, deprecatedCID_cid, deprecatedRID_rid)
        rpCache._gen_comp_xref_deprecatedCompID_compid(input_dir, outdir)
        rpCache._gen_rr_full_reactions(input_dir, outdir, deprecatedCID_cid, deprecatedRID_rid)
        del deprecatedCID_cid, deprecatedRID_rid


    @staticmethod
    def _gen_deprecatedCID_cid(input_dir, outdir):
        attribute = 'deprecatedCID_cid'
        print(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedCID_cid = None
        f_deprecatedCID_cid = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_deprecatedCID_cid):
            print("   Generating data...", end = '', flush=True)
            deprecatedCID_cid = rpCache._m_deprecatedMNXM(input_dir+'chem_xref.tsv.gz')
            #overwrite (or not if it dosn't exist) entries that are defined by Thomas
            try:
                user_mnx_replace = json_load(open('data/mnx_replace.json', 'r'))
                for user_deprecated_mnxm in user_mnx_replace:
                    deprecatedCID_cid[user_deprecated_mnxm] = user_mnx_replace[user_deprecated_mnxm]['mnx']
            except FileNotFoundError:
                print("   Error data/mnx_replace.json file not found")
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(deprecatedCID_cid, f_deprecatedCID_cid)
            print_OK()
        else:
            deprecatedCID_cid = rpCache._load_cache_from_file(f_deprecatedCID_cid)
            print("   Cache file already exists", end = '', flush=True)
            print_OK()
        return {'attr': deprecatedCID_cid, 'file': f_deprecatedCID_cid}


    @staticmethod
    def _gen_cid_strc_cid_name(input_dir, outdir, deprecatedCID_cid):
        attribute = 'cid_strc, cid_name'
        print(c_attr('bold')+attribute+c_attr('reset'))
        cid_strc = None
        cid_name = None
        f_cid_strc = outdir+'cid_strc'+rpCache._ext
        f_cid_name = outdir+'cid_name'+rpCache._ext
        if not os_path.isfile(f_cid_strc):
            if not deprecatedCID_cid['attr']:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedCID_cid = rpCache._load_cache_from_file(deprecatedCID_cid['file'])
                print_OK()
            print("   Generating data...", end = '', flush=True)
            cid_strc, cid_name = rpCache._m_mnxm_strc(input_dir+'/compounds.tsv.gz', input_dir+'chem_prop.tsv.gz', deprecatedCID_cid['attr'])
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(cid_strc, f_cid_strc)
            rpCache._store_cache_to_file(cid_name, f_cid_name)
            print_OK()
        else:
            cid_strc = rpCache._load_cache_from_file(f_cid_strc)
            print("   Cache file already exists", end = '', flush=True)
            print_OK()
        return {'attr': cid_strc, 'file': f_cid_strc}, {'attr': cid_name, 'file': f_cid_name}


    @staticmethod
    def _gen_inchikey_cid(input_dir, outdir, cid_strc):
        attribute = 'inchikey_cid'
        print(c_attr('bold')+attribute+c_attr('reset'))
        inchikey_cid = None
        f_inchikey_cid = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_inchikey_cid):
            if not cid_strc['attr']:
                print("   Loading input data from file...", end = '', flush=True)
                cid_strc['attr'] = rpCache._load_cache_from_file(cid_strc['file'])
                print_OK()
            print("   Generating data...", end = '', flush=True)
            inchikey_cid = rpCache._m_inchikey_cid(cid_strc['attr'])
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(inchikey_cid, f_inchikey_cid)
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()


    @staticmethod
    def _gen_cid_xref(input_dir, outdir, deprecatedCID_cid):
        attribute = 'cid_xref'
        print(c_attr('bold')+attribute+c_attr('reset'))
        cid_xref = None
        f_cid_xref = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_cid_xref):
            if not deprecatedCID_cid['attr']:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedCID_cid['attr'] = rpCache._load_cache_from_file(deprecatedCID_cid['file'])
                print_OK()
            print("   Generating data...", end = '', flush=True)
            cid_xref = rpCache._m_mnxm_xref(input_dir+'chem_xref.tsv.gz', deprecatedCID_cid['attr'])
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(cid_xref, f_cid_xref)
            print_OK()
        else:
            cid_xref = rpCache._load_cache_from_file(f_cid_xref)
            print("   Cache file already exists", end = '', flush=True)
            print_OK()
        return {'attr': cid_xref, 'file': f_cid_xref}


    @staticmethod
    def _gen_chebi_cid(input_dir, outdir, cid_xref):
        attribute = 'chebi_cid'
        print(c_attr('bold')+attribute+c_attr('reset'))
        chebi_cid = None
        f_chebi_cid = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_chebi_cid):
            print("   Generating data...", end = '', flush=True)
            chebi_cid = rpCache._m_chebi_cid(cid_xref['attr'])
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(chebi_cid, f_chebi_cid)
            del chebi_cid
            print_OK()
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()


    @staticmethod
    def _gen_deprecatedRID_rid(input_dir, outdir):
        attribute = 'deprecatedRID_rid'
        print(c_attr('bold')+attribute+c_attr('reset'))
        deprecatedRID_rid = None
        f_deprecatedRID_rid = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_deprecatedRID_rid):
            print("   Generating data...", end = '', flush=True)
            deprecatedRID_rid = rpCache._m_deprecatedMNXR(input_dir+'reac_xref.tsv.gz')
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(deprecatedRID_rid, f_deprecatedRID_rid)
            print_OK()
        else:
            deprecatedRID_rid = rpCache._load_cache_from_file(f_deprecatedRID_rid)
            print("   Cache file already exists", end = '', flush=True)
            print_OK()
        return {'attr': deprecatedRID_rid, 'file': f_deprecatedRID_rid}


    @staticmethod
    def _gen_rr_reactions(input_dir, outdir, deprecatedCID_cid, deprecatedRID_rid):
        attribute = 'rr_reactions'
        print(c_attr('bold')+attribute+c_attr('reset'))
        rr_reactions = None
        f_rr_reactions = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_rr_reactions):
            if not deprecatedCID_cid['attr']:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedCID_cid['attr'] = rpCache._load_cache_from_file(deprecatedCID_cid['file'])
                print_OK()
            if not deprecatedRID_rid['attr']:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedRID_rid['attr'] = rpCache._load_cache_from_file(deprecatedRID_rid['file'])
                print_OK()
            print("   Generating data...", end = '', flush=True)
            rr_reactions = rpCache._m_rr_reactions(input_dir+'retrorules_rr02_flat_all.tsv.gz', deprecatedCID_cid, deprecatedRID_rid)
            print_OK()
            del deprecatedRID_rid
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(rr_reactions, f_rr_reactions)
            print_OK()
            del rr_reactions
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()
        # return deprecatedCID_cid


    @staticmethod
    def _gen_comp_xref_deprecatedCompID_compid(input_dir, outdir):
        attribute = 'comp_xref, deprecatedCompID_compid'
        print(c_attr('bold')+attribute+c_attr('reset'))
        comp_xref = deprecatedCompID_compid = None
        f_comp_xref = outdir+'comp_xref'+rpCache._ext
        f_deprecatedCompID_compid = outdir+'deprecatedCompID_compid'+rpCache._ext
        if not os_path.isfile(f_comp_xref) or not os_path.isfile(f_deprecatedCompID_compid):
            print("   Generating data...", end = '', flush=True)
            comp_xref,deprecatedCompID_compid = rpCache._m_mnxc_xref(input_dir+'comp_xref.tsv.gz')
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(comp_xref, f_comp_xref)
            print_OK()
            del comp_xref
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(deprecatedCompID_compid, f_deprecatedCompID_compid)
            print_OK()
            del deprecatedCompID_compid
        else:
            print("   Cache files already exist", end = '', flush=True)
            print_OK()


    @staticmethod
    def _gen_rr_full_reactions(input_dir, outdir, deprecatedCID_cid, deprecatedRID_rid):
        attribute = 'rr_full_reactions'
        print(c_attr('bold')+attribute+c_attr('reset'))
        rr_full_reactions = None
        f_rr_full_reactions = outdir+attribute+rpCache._ext
        if not os_path.isfile(f_rr_full_reactions):
            print("   Generating data...", end = '', flush=True)
            if not deprecatedCID_cid['attr']:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedCID_cid = rpCache._load_cache_from_file(deprecatedCID_cid['file'])
                print_OK()
            if not deprecatedRID_rid:
                print("   Loading input data from file...", end = '', flush=True)
                deprecatedRID_rid = rpCache._load_cache_from_file(deprecatedRID_rid['file'])
                print_OK()
            rr_full_reactions = rpCache._m_rr_full_reactions(input_dir+'rxn_recipes.tsv.gz', deprecatedCID_cid['attr'], deprecatedRID_rid['attr'])
            print_OK()
            print("   Writing data to file...", end = '', flush=True)
            rpCache._store_cache_to_file(rr_full_reactions, f_rr_full_reactions)
            print_OK()
            del rr_full_reactions
        else:
            print("   Cache file already exists", end = '', flush=True)
            print_OK()


    @staticmethod
    def _check_or_download_cache_to_disk(cache_dir, attributes):
        for attr in attributes:
            filename = attr+rpCache._ext
            if os_path.isfile(cache_dir+filename) and sha512(Path(cache_dir+filename).read_bytes()).hexdigest()==rpCache._cache_files[filename]:
                print(filename+" already downloaded ", end = '', flush=True)
                print_OK()
            else:
                filename = attr+rpCache._ext
                print("Downloading "+filename+"...", end = '', flush=True)
                start_time = time_time()
                if not os_path.isdir(cache_dir):
                    os_mkdir(cache_dir)
                download(rpCache._cache_url+filename, cache_dir+filename)
                rpCache._cache_files[attr] = True
                end_time = time_time()
                print_OK(end_time-start_time)


    def _load_from_file(self, attribute):
        filename = attribute+rpCache._ext
        print("Loading "+filename+"...", end = '', flush=True)
        data = self._load_cache_from_file(self._cache_dir+filename)
        print_OK()
        return data


    def _check_or_load_cache(self):
        if self.store_mode=='file':
            self._check_or_load_cache_in_memory()
        else:
            self._check_or_load_cache_in_db()


    def _check_or_load_cache_in_memory(self):
        for attribute in self._attributes:
            if not getattr(self, attribute):
                setattr(self, attribute, self._load_from_file(attribute))
            else:
                print(attribute+" already loaded in memory...", end = '', flush=True)
                print_OK()


    def _check_or_load_cache_in_db(self):
        for attribute in self._attributes:
            if not CRedisDict.exists(self.redis, attribute):
                self._store_cache_to_db(attribute, self._load_from_file(attribute))
            else:
                print(attribute+" already loaded in db...", end = '', flush=True)
                print_OK()


    @staticmethod
    def _download_input_cache(url, file, outdir):
        if not os_path.isdir(outdir):
            os_mkdir(outdir)
        filename = outdir+file
        if not os_path.isfile(filename):
            print("Downloading "+file+"...", end = '', flush=True)
            start_time = time_time()
            rpCache.__download_input_cache(url, file, outdir)
            end_time = time_time()
            print_OK(end_time-start_time)
        else:
            print(filename+" already downloaded ", end = '', flush=True)
            print_OK()


    @staticmethod
    def __download_input_cache(url, file, outdir):

        if not os_path.isdir(outdir):
            os_mkdir(outdir)


        # 3xCommon + rpReader
        if file in ['reac_xref.tsv.gz', 'chem_xref.tsv.gz', 'chem_prop.tsv.gz', 'comp_xref.tsv.gz']:
            download(url+'metanetx/'+file, outdir+file)

        #TODO: need to add this file to the git or another location
        if file in ['compounds.tsv.gz', 'rxn_recipes.tsv.gz']:
            download(url+'rr02_more_data/'+file,
                     outdir+file)
            # tar = tarfile_open(outdir+'/rr02_more_data.tar.gz', 'r:gz')
            # tar.extractall(outdir)
            # tar.close()
            # shutil_move(outdir+'/rr02_more_data/compounds.tsv',
            #             outdir+'/rr_compounds.tsv')
            # shutil_move(outdir+'/rr02_more_data/rxn_recipes.tsv',
            #             outdir)
            # os_rm(outdir+'/rr02_more_data.tar.gz')
            # shutil_rmtree(outdir+'/rr02_more_data')

        if file=='retrorules_rr02_flat_all.tsv.gz':
            download(url+'retrorules_rr02_rp3_hs/'+file,
                     outdir+file)
            # download('https://retrorules.org/dl/preparsed/rr02/rp3/hs',
            #          outdir+'/retrorules_rr02_rp3_hs.tar.gz')
            # tar = tarfile_open(outdir+'/retrorules_rr02_rp3_hs.tar.gz', 'r:gz')
            # tar.extractall(outdir)
            # tar.close()
            # shutil_move(outdir+'/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv', outdir+'/rules_rall.tsv')
            # os_rm(outdir+'/retrorules_rr02_rp3_hs.tar.gz')
            # shutil_rmtree(outdir+'/retrorules_rr02_rp3_hs')

    ##########################################################
    ################## Private Functions #####################
    ##########################################################

    ## Method to load data from file
    #
    #  Load data from file
    #
    #  @param self Object pointer
    #  @param filename File to fetch data from
    #  @return file content
    @staticmethod
    def _load_cache_from_file(filename):
        if filename.endswith('.gz') or filename.endswith('.zip'):
            fp = gzip_open(filename, 'rt', encoding='ascii')
        else:
            fp = open(filename, 'r')
        return json_load(fp)

    ## Method to store data into file
    #
    # Store data into file as json (to store dictionnary structure)
    #
    #  @param self Object pointer
    #  @param data Data to write into file
    #  @param filename File to write data into
    @staticmethod
    def _store_cache_to_file(data, filename):
        if filename.endswith('.gz') or filename.endswith('.zip'):
            fp = gzip_open(filename, 'wt', encoding='ascii')
        else:
            fp = open(filename, 'w')
        json_dump(data, fp)

    ## Method to store data into redis database
    #
    #  Assign a CRedisDict object to the attribute to copy data into the database
    #
    #  @param self Object pointer
    #  @param attr_name Attribute name (database key)
    #  @param data Content of the attribute
    def _store_cache_to_db(self, attr_name, data):
        print("Storing "+attr_name+" to db...", end = '', flush=True)
        setattr(rpCache, attr_name, CRedisDict(attr_name, self.redis, data))
        print_OK()



    ## Function to create a dictionnary of old to new chemical id's
    #
    #  Generate a one-to-one dictionnary of old id's to new ones. Private function
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    @staticmethod
    def _checkCIDdeprecated(cid, deprecatedCID_cid):
        try:
            return deprecatedCID_cid[cid]
        except (KeyError, TypeError):
            return cid


    ## Function to create a dictionnary of old to new reaction id's
    #
    # TODO: check other things about the mnxm emtry like if it has the right structure etc...
    @staticmethod
    def _checkRIDdeprecated(rid, deprecatedRID_rid):
        try:
            return deprecatedRID_rid[rid]
        except (KeyError, TypeError):
            return rid


    #################################################################
    ################## Public functions #############################
    #################################################################


    ########################### MNX parsers #############################

    ## Function to parse the chem_xref.tsv and reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    # @param xref_path Input file path
    # @return Dictionnary of identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    @staticmethod
    def _deprecatedMNX(xref_path):
        deprecatedMNX_mnx = {}
        with gzip_open(xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = row[0].split(':')
                    if mnx[0]=='deprecated':
                        deprecatedMNX_mnx[mnx[1]] = row[1]
        return deprecatedMNX_mnx

    ## Status function that parses the chem_xref.tsv file for chemical cross-references and the different
    #
    # @param chem_xref_path Input file path
    # @return Dictionnary of chemical id to other chemical ids ex: deprecatedCID_cid['MNXM1'] = {'mnx': ['MNXM01', ...], ...}
    @staticmethod
    def _m_deprecatedMNXM(chem_xref_path):
        deprecatedCID_cid = {}
        deprecatedCID_cid = rpCache._deprecatedMNX(chem_xref_path)
        deprecatedCID_cid.update(rpCache._convertMNXM)
        deprecatedCID_cid['MNXM01'] = 'MNXM1'
        return deprecatedCID_cid

    ## Function to parse the reac_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of old to new MetanetX identifiers to make sure that we always use the freshest id's.
    # This can include more than one old id per new one and thus returns a dictionnary. Private function
    #
    #  @param self Object pointer
    #  @param reac_xref_path Input file path
    #  @return Dictionnary of identifiers
    @staticmethod
    def _m_deprecatedMNXR(reac_xref_path):
        return rpCache._deprecatedMNX(reac_xref_path)


    ## Function to parse the chemp_prop.tsv file from MetanetX and compounds.tsv from RetroRules. Uses the InchIkey as key to the dictionnary
    #
    #  Generate a dictionnary gaving the formula, smiles, inchi and inchikey for the components
    # TODO: Seperate this function to parse the chem_prop (mnx specific) and the compounds.tsv from RetroRules (generic, not mnx specific)
    # Structure of return: cid_strc['MNXM1'] = {'formula': 'H', 'smiles': '[H+]', 'inchi': 'InChI=1S/p+1', 'inchikey': 'GPRLSGONYQIRFK-UHFFFAOYSA-N'}
    #
    #  @param rr_compounds_path Path to the RetroRules file
    #  @param chem_prop_path Path to the chem_prop.tsv file
    #  @param deprecatedCID_cid Dictionnary of deprecated CID to cid
    #  @return cid_strc Dictionnary of formula, smiles, inchi and inchikey
    @staticmethod
    def _m_mnxm_strc(rr_compounds_path, chem_prop_path, deprecatedCID_cid, logger=None):
        logger = logger or logging.getLogger(__name__)
        cid_strc = {}
        cid_name = {}
        for row in csv_DictReader(gzip_open(rr_compounds_path, 'rt'), delimiter='\t'):
            tmp = {'formula':  None,
                    'smiles': None,
                    'inchi': row['inchi'],
                    'inchikey': None,
                    'cid': rpCache._checkCIDdeprecated(row['cid'], deprecatedCID_cid),
                    'name': None}
            try:
                resConv = rpCache._convert_depiction(idepic=tmp['inchi'], itype='inchi', otype={'smiles','inchikey'})
                for i in resConv:
                    tmp[i] = resConv[i]
            except rpCache.DepictionError as e:
                logger.warning('Could not convert some of the structures: '+str(tmp))
                logger.warning(e)
            cid_strc[tmp['cid']] = tmp
        with gzip_open(chem_prop_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnxm = rpCache._checkCIDdeprecated(row[0], deprecatedCID_cid)
                    tmp = {'formula':  row[2],
                            'smiles': row[6],
                            'inchi': row[5],
                            'inchikey': row[8],
                            'cid': mnxm,
                            'name': row[1]}
                    for i in tmp:
                        if tmp[i]=='' or tmp[i]=='NA':
                            tmp[i] = None
                    if not mnxm in cid_name and tmp['name']:
                        cid_name[mnxm] = tmp['name']
                    if mnxm in cid_strc:
                        cid_strc[mnxm]['formula'] = row[2]
                        cid_strc[mnxm]['name'] = row[1]
                        if not cid_strc[mnxm]['smiles'] and tmp['smiles']:
                            cid_strc[mnxm]['smiles'] = tmp['smiles']
                        if not cid_strc[mnxm]['inchikey'] and tmp['inchikey']:
                            cid_strc[mnxm]['inchikey'] = tmp['inchikey']
                    else:
                        #check to see if the inchikey is valid or not
                        otype = set({})
                        if not tmp['inchikey']:
                            otype.add('inchikey')
                        if not tmp['smiles']:
                            otype.add('smiles')
                        if not tmp['inchi']:
                            otype.add('inchi')
                        itype = ''
                        if tmp['inchi']:
                            itype = 'inchi'
                        elif tmp['smiles']:
                            itype = 'smiles'
                        else:
                            logger.warning('No valid entry for the convert_depiction function')
                            continue
                        try:
                            resConv = rpCache._convert_depiction(idepic=tmp[itype], itype=itype, otype=otype)
                            for i in resConv:
                                tmp[i] = resConv[i]
                        except rpCache.DepictionError as e:
                            logger.warning('Could not convert some of the structures: '+str(tmp))
                            logger.warning(e)
                        cid_strc[tmp['cid']] = tmp
        return cid_strc, cid_name


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #
    #  @param chem_xref_path MetaNetX chem_xref.tsv file path
    #  @param deprecatedCID_cid Dictionnary of deprecated chemical ids to uniform cid
    #  @return Dictionnary of cross references of a given chemical id
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    @staticmethod
    def _m_mnxm_xref(chem_xref_path, deprecatedCID_cid, logger=None):
        logger = logger or logging.getLogger(__name__)
        cid_xref = {}
        with gzip_open(chem_xref_path, 'rt') as f:
            c = csv_reader(f, delimiter='\t')
            for row in c:
                if not row[0][0]=='#':
                    mnx = rpCache._checkCIDdeprecated(row[1], deprecatedCID_cid)
                    if len(row[0].split(':'))==1:
                        dbName = 'mnx'
                        dbId = row[0]
                    else:
                        dbName = row[0].split(':')[0]
                        dbId = ''.join(row[0].split(':')[1:])
                        if dbName=='deprecated':
                            dbName = 'mnx'
                    #mnx
                    if not mnx in cid_xref:
                        cid_xref[mnx] = {}
                    if not dbName in cid_xref[mnx]:
                        cid_xref[mnx][dbName] = []
                    if not dbId in cid_xref[mnx][dbName]:
                        cid_xref[mnx][dbName].append(dbId)
                    ### DB ###
                    if not dbName in cid_xref:
                        cid_xref[dbName] = {}
                    if not dbId in cid_xref[dbName]:
                        cid_xref[dbName][dbId] = mnx
        return cid_xref


    ## Function to parse the comp_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of compartments id's (MNX) to other database id's
    #
    #  @param comp_xref_path The MetaNetX file that contains the cross references
    #  @return a The dictionnary of compartment identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
    @staticmethod
    def _m_mnxc_xref(comp_xref_path, logger=None):
        logger = logger or logging.getLogger(__name__)
        comp_xref = {}
        deprecatedCompID_compid = {}
        try:
            with gzip_open(comp_xref_path, 'rt') as f:
                c = csv_reader(f, delimiter='\t')
                #not_recognised = []
                for row in c:
                    #cid = row[0].split(':')
                    if not row[0][0]=='#':
                        #collect the info
                        mnxc = row[1]
                        if len(row[0].split(':'))==1:
                            dbName = 'mnx'
                            dbCompId = row[0]
                        else:
                            dbName = row[0].split(':')[0]
                            dbCompId = ''.join(row[0].split(':')[1:])
                            dbCompId = dbCompId.lower()
                        if dbName=='deprecated':
                            dbName = 'mnx'
                        #create the dicts
                        if not mnxc in comp_xref:
                            comp_xref[mnxc] = {}
                        if not dbName in comp_xref[mnxc]:
                            comp_xref[mnxc][dbName] = []
                        if not dbCompId in comp_xref[mnxc][dbName]:
                            comp_xref[mnxc][dbName].append(dbCompId)
                        #create the reverse dict
                        if not dbCompId in deprecatedCompID_compid:
                            deprecatedCompID_compid[dbCompId] = mnxc
        except FileNotFoundError:
            logger.error('comp_xref file not found')
            return {}
        return comp_xref,deprecatedCompID_compid


    ######################## RetroRules specific functions ##################


    ## Function to parse the rules_rall.tsv from RetroRules
    #
    #  Extract from the reactions rules the ruleID, the reactionID, the direction of the rule directed to the origin reaction
    #  Structure of the return: rr_reactions['RR-02-d2e7c5761b5a9b4b-04-F'] = {'MNXR139133': {'rule_id': 'RR-02-d2e7c5761b5a9b4b-04-F', 'rule_score': 0.3151075983206353, 'reac_id': 'MNXR139133', 'subs_id': 'MNXM89557', 'rel_direction': 1, 'left': {'MNXM89557': 1}, 'right': {'MNXM20': 1, 'MNXM722724': 1}}}
    #
    #  @param rules_rall_path Path to the RetroRules reaction rules
    #  @param deprecatedCID_cid Dictionnary of deprecated to uniformed chemical id's
    #  @param deprecatedRID_rid Dictionnary of deprecated to uniformed reaction id's
    #  @return Dictionnary describing each reaction rule
    @staticmethod
    def _m_rr_reactions(rules_rall_path, deprecatedCID_cid, deprecatedRID_rid, logger=None):
        logger = logger or logging.getLogger(__name__)
        rr_reactions = {}
        try:
            #with gzip_open(rules_rall_path, 'r') as f:
            #    reader = csv.reader(f, delimiter = '\t')
            #    next(reader)
            #    rule = {}
            #    for row in reader:
            for row in csv_DictReader(gzip_open(rules_rall_path, 'rt'), delimiter='\t'):
                #NOTE: as of now all the rules are generated using MNX
                #but it may be that other db are used, we are handling this case
                #WARNING: can have multiple products so need to seperate them
                products = {}
                for i in row['Product_IDs'].split('.'):
                    cid = rpCache._checkCIDdeprecated(i, deprecatedCID_cid)
                    if not cid in products:
                        products[cid] = 1
                    else:
                        products[cid] += 1
                try:
                    #WARNING: one reaction rule can have multiple reactions associated with them
                    #To change when you can set subpaths from the mutliple numbers of
                    #we assume that the reaction rule has multiple unique reactions associated
                    if row['# Rule_ID'] not in rr_reactions:
                        rr_reactions[row['# Rule_ID']] = {}
                    if row['# Rule_ID'] in rr_reactions[row['# Rule_ID']]:
                        logger.warning('There is already reaction '+str(row['# Rule_ID'])+' in reaction rule '+str(row['# Rule_ID']))
                    rr_reactions[row['# Rule_ID']][row['Reaction_ID']] = {
                        'rule_id': row['# Rule_ID'],
                        'rule_score': float(row['Score_normalized']),
                        'reac_id': rpCache._checkRIDdeprecated(row['Reaction_ID'], deprecatedRID_rid),
                        'subs_id': rpCache._checkCIDdeprecated(row['Substrate_ID'], deprecatedCID_cid),
                        'rel_direction': int(row['Rule_relative_direction']),
                        'left': {rpCache._checkCIDdeprecated(row['Substrate_ID'], deprecatedCID_cid): 1},
                        'right': products}
                except ValueError:
                    logger.error('Problem converting rel_direction: '+str(row['Rule_relative_direction']))
                    logger.error('Problem converting rule_score: '+str(row['Score_normalized']))
            return rr_reactions
        except FileNotFoundError as e:
                logger.error('Could not read the rules_rall file ('+str(rules_rall_path)+')')
                return {}


    ## Generate complete reactions from the rxn_recipes.tsv from RetroRules
    #
    #  These are the compplete reactions from which the reaction rules are generated from. This is used to
    #  reconstruct the full reactions from monocomponent reactions
    #  Structur of the return: rr_full_reactions['MNXR142257'] = {'left': {'MNXM4660': 1}, 'right': {'MNXM97172': 1}, 'direction': 0, 'main_left': ['MNXM4660'], 'main_right': ['MNXM97172']}
    #
    #  @param self The pointer object
    #  @param rxn_recipes_path Path to the recipes file
    #  @return Boolean that determines the success or failure of the function
    @staticmethod
    def _m_rr_full_reactions(rxn_recipes_path, deprecatedCID_cid, deprecatedRID_rid, logger=None):
        logger = logger or logging.getLogger(__name__)
        #### for character matching that are returned
        DEFAULT_STOICHIO_RESCUE = {"4n": 4, "3n": 3, "2n": 2, 'n': 1,
                           '(n)': 1, '(N)': 1, '(2n)': 2, '(x)': 1,
                           'N': 1, 'm': 1, 'q': 1,
                           '0.01': 1, '0.1': 1, '0.5': 1, '1.5': 1,
                           '0.02': 1, '0.2': 1,
                           '(n-1)': 0, '(n-2)': -1}
        reaction = {}
        try:
            for row in csv_DictReader(gzip_open(rxn_recipes_path, 'rt'), delimiter='\t'):
                tmp = {} # makes sure that if theres an error its not added
                #parse the reaction equation
                if not len(row['Equation'].split('='))==2:
                    logger.warning('There should never be more or less than a left and right of an equation')
                    logger.warnin(row['Equation'])
                    continue
                ######### LEFT ######
                #### MNX id
                tmp['left'] = {}
                # if row['#Reaction_ID']=="MNXR141948":
                #     print(row)
                #     exit()
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[0]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['left'][rpCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['left'][rpCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = int(spe[0])
                        except ValueError:
                            logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### RIGHT #####
                ####  MNX id
                tmp['right'] = {}
                for spe in re_findall(r'(\(n-1\)|\d+|4n|3n|2n|n|\(n\)|\(N\)|\(2n\)|\(x\)|N|m|q|\(n\-2\)|\d+\.\d+) ([\w\d]+)@\w+', row['Equation'].split('=')[1]):
                    #1) try to rescue if its one of the values
                    try:
                        tmp['right'][rpCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = DEFAULT_STOICHIO_RESCUE[spe[0]]
                    except KeyError:
                        #2) try to convert to int if its not
                        try:
                            tmp['right'][rpCache._checkCIDdeprecated(spe[1], deprecatedCID_cid)] = int(spe[0])
                        except ValueError:
                            logger.warning('Cannot convert '+str(spe[0]))
                            continue
                ####### DIRECTION ######
                try:
                    tmp['direction'] = int(row['Direction'])
                except ValueError:
                    logger.error('Cannot convert '+str(row['Direction'])+' to int')
                    continue
                ### add the others
                tmp['main_left'] = row['Main_left'].split(',')
                tmp['main_right'] = row['Main_right'].split(',')
                reaction[rpCache._checkRIDdeprecated(row['#Reaction_ID'], deprecatedRID_rid)] = tmp
            return reaction
        except FileNotFoundError:
            logger.error('Cannot find file: '+str(rxn_recipes_path))
            return False


    ######################## Generic functions ###############################

    ## Convert chemical depiction to others type of depictions
    #
    # Usage example:
    # - convert_depiction(idepic='CCO', otype={'inchi', 'smiles', 'inchikey'})
    # - convert_depiction(idepic='InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', itype='inchi', otype={'inchi', 'smiles', 'inchikey'})
    #
    #  @param self The object pointer
    #  @param idepic String depiction to be converted, str
    #  @param itype type of depiction provided as input, str
    #  @param otype types of depiction to be generated, {"", "", ..}
    #  @return odepic generated depictions, {"otype1": "odepic1", ..}
    @staticmethod
    def _convert_depiction(idepic, itype='smiles', otype={'inchikey'}):
        # Import (if needed)
        if itype == 'smiles':
            rdmol = MolFromSmiles(idepic, sanitize=True)
        elif itype == 'inchi':
            rdmol = MolFromInchi(idepic, sanitize=True)
        else:
            raise NotImplementedError('"{}" is not a valid input type'.format(itype))
        if rdmol is None:  # Check imprt
            raise rpCache.DepictionError('Import error from depiction "{}" of type "{}"'.format(idepic, itype))
        # Export
        odepic = dict()
        for item in otype:
            if item == 'smiles':
                odepic[item] = MolToSmiles(rdmol)  # MolToSmiles is tricky, one mays want to check the possible options..
            elif item == 'inchi':
                odepic[item] = MolToInchi(rdmol)
            elif item == 'inchikey':
                odepic[item] = MolToInchiKey(rdmol)
            else:
                raise NotImplementedError('"{}" is not a valid output type'.format(otype))
        return odepic


    ## Function to parse the chem_xref.tsv file of MetanetX
    #
    #  Generate a dictionnary of all cross references for a given chemical id (MNX) to other database id's
    #  Structure if the return: chebi_cid['88281']: 'MXM2323'
    #
    #  @param self Object pointer
    #  @param chem_xref_path Input file path
    #  @return a The dictionnary of identifiers
    #TODO: save the self.deprecatedCID_cid to be used in case there rp_paths uses an old version of MNX
#    def _m_chebi_cid(self, cid_xref):
    @staticmethod
    def _m_chebi_cid(cid_xref):
        chebi_cid = {}
        for cid in cid_xref:
            if 'chebi' in cid_xref[cid]:
                for c in cid_xref[cid]['chebi']:
                    chebi_cid[c] = cid
        return chebi_cid

    ## Function to build the dictionnary to find the chemical id from inchikey
    #
    # @param cid_strc Dictionnary of chemical ID's to all the structure information associated with it
    # @return Dictionnary of InChIKey to chemical ID
    @staticmethod
    def _m_inchikey_cid(cid_strc):
        inchikey_cid = {}
        for cid in cid_strc:
            inchikey = cid_strc[cid]['inchikey']
            # This line is needed to put a value in 'inchikey', otherwise there are some problems in future strucutres
            if not inchikey: inchikey = 'NO_INCHIKEY'
            if not inchikey in inchikey_cid:
                inchikey_cid[inchikey] = []
            inchikey_cid[inchikey].append(cid)
        return inchikey_cid
