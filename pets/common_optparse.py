class Common_optparse:
    def add_options(parser):
        parser.add_argument("-c", "--chromosome_col", dest="chromosome_col", default= None,
                            help="Column name if header is true, otherwise 0-based position of the column with the chromosome")

        parser.add_argument("-d", "--pat_id_col", dest="id_col", default= None,
                            help="Column name if header is true, otherwise 0-based position of the column with the patient id")
        
        parser.add_argument("-s", "--start_col", dest="start_col", default= None,
                            help="Column name if header is true, otherwise 0-based position of the column with the start mutation coordinate")
        
        parser.add_argument("-e", "--end_col", dest="end_col", default= None,
                            help="Column name if header is true, otherwise 0-based position of the column with the end mutation coordinate")

        parser.add_argument("-G", "--genome_assembly", dest="genome_assembly", default= "hg38",
                            help="Genome assembly version. Please choose between hg18, hg19 and hg38. Default hg38")

        #chr\tstart\tstop
        parser.add_argument("-H", "--header", dest="header", default= True, action="store_false",
                            help="Set if the file has a line header. Default true")

        parser.add_argument("-x", "--sex_col", dest="sex_col", default= None,
                            help="Column name if header is true, otherwise 0-based position of the column with the patient sex")

        parser.add_argument("-p", "--hpo_term_col", dest="ont_col", default= None,
                            help="Column name if header true or 0-based position of the column with the HPO terms")
        
        parser.add_argument("-S", "--hpo_separator", dest="separator", default='|',
                            help="Set which character must be used to split the HPO profile. Default '|'")

        parser.add_argument("-n","--hpo_names", dest="names", default= False, action="store_true",
                            help="Define if the input HPO are human readable names. Default false")

        parser.add_argument("-X", "--excluded_hpo", dest="excluded_hpo", default= None,
                            help="File with excluded HPO terms")