import argparse as ap

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")  
    parser.add_argument('-sample', required=True, help='sample list for cohort')
    parser.add_argument('-bim', required=True, help='variant file')
    parser.add_argument('-per_gpu_mem', required=True, help='max memory per gpu', default=10)
    parser.add_argument('-c', '--cohort', required=True, help='Cohort')
    parser.add_argument('-p', '--pheno', help='Phenotype')

    return parser


args = make_arg_parser().parse_args()
n_var = args.bim
n_sample = args.sample

approx_mem_usage = 4 * n_var * n_sample / 1e9
n_gpu_per = int(np.ceil(approx_mem_usage / args.per_gpu_mem))