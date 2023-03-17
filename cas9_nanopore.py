#!/usr/bin/env python
# coding: utf-8
#主程序源代码部分，配置文件部分在文末
import os
import sys
import shutil
import argparse
import configparser
import logging
import time
def filter(fastq , outFq , minlength ):
    #按reads长度多虑fastq文件
    with open(fastq,"r") as fqfile ,open(outFq,'w') as fqfilter:
        i=0
        for line in fqfile:
            if i%4==0:
                seqID=line
            elif i%4==1:
                sequence=line
            elif i%4==2:
                mark=line
            elif i%4==3:
                if len(sequence)>=minlength:
                    qual=line.strip("\n")
                    fqfilter.write(seqID+sequence+mark+qual+"\n")
            i+=1
class megalodon(object):
    def __init__(self ,parser,argv, config):
        self.argv = argv
        self.mean = ''
        self.cmd = ''
        self.args = ''
        self.parser = parser
        self.config = config
    def extract_param(self):
        #设置输入参数，部分参数在配置文件更改
        parser = self.parser
        parser.add_argument('-i', '--input', type=str, required=True, help='input the FAST5 folder')
        parser.add_argument('-o', '--output', type=str, required=True, help='output folder')
        parser.add_argument('-g', '--genome', type=str, required = True, help='reference genome sequences in fasta format')
        parser.add_argument('--guppy_config', default='',
                            help=' Guppy config.'
                                 'Change the path in the config.ini file.')
        parser.add_argument('--guppy_server_path',type=str, default = '',
                            help='Absolute path to the guppy_basecall_server.\
                                 The default value is in the config.ini.')
        parser.add_argument('-b','--trim_barcode', action='store_true', default=False,
                            help='Whether to remove the barcode sequence.\
                                 default = False')
        parser.add_argument('--remora', action='store_true', default=False,
                            help='Whether to use the remora model.\
                                 default: rerio model')
        parser.add_argument('-t', type=int, default=16,
                            help='Number of parallel processes.default:16')
        parser.add_argument('--guppy_param', default='',
                            help='The default value is in the config.ini.\
                                 This parameter cannot be used in remora mode')
        self.args =  parser.parse_args(self.argv)
    def config_praram(self):
        #从配置文件中获取部分参数
        args = self.args
        if args.output[-1] == '/':
            args.output = args.output[:-1]
        config = self.config
        if args.remora:
            if not args.guppy_config:
                args.guppy_config = config.get('remora', 'guppy_config')
            if not args.guppy_server_path:
                args.guppy_server_path = config.get('remora','guppy_server_path')
        else:
            if not args.guppy_param:
                args.guppy_param = config.get('rerio', 'guppy_param')
            if not args.guppy_config:
                args.guppy_config = config.get('rerio', 'guppy_config')
            if not args.guppy_server_path:
                args.guppy_server_path = config.get('rerio','guppy_server_path')
        if args.trim_barcode:
            args.guppy_param += ' --trim_barcode'
    def run_cmd(self):
        args = self.args
        if args.remora:
            self.means = 'remora'
            self.cmd = 'megalodon {input_fast5} \
--guppy-config {cfg} \
--remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmC_5mC CG 0 \
--outputs basecalls mappings mod_mappings mods per_read_mods \
--reference {ref} \
--devices 0 \
--processes {t} \
--mod-map-emulate-bisulfite --mod-map-base-conv C T \
--mod-map-base-conv m C \
--guppy-server-path {path} \
--output-directory {output} \
--write-mods-text' \
                .format(input_fast5=args.input,
                        cfg=args.guppy_config,
                        ref=args.genome,
                        t=args.t,
                        path=args.guppy_server_path,
                        output=args.output)
        else:
            self.means = 'rerio'
            self.cmd = 'megalodon {input_fast5} --guppy-param \'{param}\' \
--guppy-config {cfg} \
--outputs basecalls mappings mod_mappings mods per_read_mods \
--reference {ref} \
--mod-motif m CG 0 \
--devices 0 \
--processes {t} \
--mod-map-emulate-bisulfite --mod-map-base-conv C T \
--mod-map-base-conv m C \
--guppy-server-path {path} \
--output-directory {output} \
--write-mods-text' \
                .format(input_fast5=args.input,
                        param=args.guppy_param,
                        cfg=args.guppy_config,
                        ref=args.genome,
                        t=args.t,
                        path=args.guppy_server_path,
                        output=args.output)
    def run_logging(self):
        args = self.args
        logging.info('\n\
megalodon analysis\n\
model:              {m}\n\
trim_barcode:       {trim}\n\
process:            {t}\n\
cmd:\n{c}'.format(m=self.means,
                          trim=args.trim_barcode,
                          t=args.t,
                          c=self.cmd))
        os.system(self.cmd)
class bcftools(object):
    def __init__(self ,parser,argv, config):
        self.argv = argv
        self.hg = ''
        self.SnpSift = ''
        self.snpEff = ''
        self.vcf_annotate = ''
        self.snpEff_config = ''
        self.args = ''
        self.parser = parser
        self.config = config
    def extract_param(self):
        parser = self.parser
        parser.add_argument('-i', '--input', required=True, help='input reads in fastq format')
        parser.add_argument('-o', '--output', required=True,
                            help="Name your project output directory.")
        parser.add_argument('-g', '--genome', required=True,
                            help="Reference genome sequence in FASTA format.\n")
        parser.add_argument('--snpEff',default='',
                            help="The path of snpEff.jar")
        parser.add_argument('--snpSift',default='',
                            help='The path of snpSift.jar')
        parser.add_argument('--snpEff_config',default='',
                            help='The path of snpEff.config')
        parser.add_argument('--snpEff_mode',default='',
                            help='Verbose mode, default:GRCh38.p14')
        parser.add_argument('--annotate',default='',
                            help='Address of the vcf annotation file')
        parser.add_argument('-t', type=int, default=16,
                            help='Number of parallel processes. default:16')
        parser.add_argument('-p', action='store_true', default=False,
                            help="Whether to use Porechop to  remove the adapter sequence."
                                 " default： False")
        parser.add_argument('-q', type=int, default=None,
                            help='minimum mapping quality of sub-reads. default:None')
        parser.add_argument('-l', '--minLength', type=int,
                            required=False, default=None, nargs='*',
                            help="Filter minimum length,Multiple values can be entered ,No filtering by default")
        parser.add_argument('-qual', type=int, default=20,
                            help='The vcf file minimum quality. default:20')
        parser.add_argument('--minDP', type=int, default=10,
                            help='The vcf file minimum depth to call variants. default：10')
        parser.add_argument('--type', type=str, default='snp',
                            help='Type of variation analyzed. default:\'snp\'')
        self.args =  parser.parse_args(self.argv)
    def config_praram(self):
        args = self.args
        config = self.config
        if args.snpEff_mode:
            self.hg = args.snpEff_mode
        else:
            self.hg = config.get('bcftools','hg')
        if args.snpSift:
            self.SnpSift = args.snpSift
        else:
            self.SnpSift = config.get('bcftools', 'SnpSift')
        if args.snpEff:
            self.snpEff = arge.snpEff
        else:
            self.snpEff = config.get('bcftools', 'snpEff')
        if args.annotate:
            self.vcf_annotate = args.annotate
        else:
            self.vcf_annotate = config.get('bcftools', 'annotate')
        if args.snpEff_config:
            self.snpEff_config = args.snpEff_config
        else:
            self.snpEff_config = config.get('bcftools', 'snpEff_config')
    def log(self):
        args = self.args
        logging.info('\n\
        input :                {i}\n\
        output :               {o}\n\
        reference :            {g}\n\
        config reaference :    {hg}\n\
        threads :              {t}\n\
        PoreChop :             {p}\n\
        minimum length :       {l}\n\
        minimum mappingQual :  {q}\n\
        minimum qual :         {qual}\n\
        minimum depth :        {dp}\n\
        variation :            {type}'.format(i=args.input,
                                              o=args.output,
                                              g=args.genome,
                                              hg=self.hg,
                                              t=args.t,
                                              p=args.p,
                                              l=args.minLength,
                                              q=args.q,
                                              qual=args.qual,
                                              dp=args.minDP,
                                              type=args.type))
    def run_cmd(self):
        args = self.args
        if args.output[-1] == '/':
            args.output = args.output[:-1]
        SnpSift = self.SnpSift
        snpEff = self.snpEff
        snpEff_config = self.snpEff_config
        vcf_annotate = self.vcf_annotate
        if not os.path.exists(args.output):
            os.mkdir(args.output)
        if args.p:
            output_porechop = args.output + '/trimmed_adapter.fastq'
            cmd = 'porechop -t {t} -i {input} -o {output}'.format(t=args.t,
                                                                  input=args.input,
                                                                  output=output_porechop)
            os.system(cmd)
            args.input = output_porechop
        if args.minLength:
            l = args.minLength
            for i in l:
                #过滤长度各输出储存一个文件夹
                if len(l) >= 1:
                    path = args.output + '/{0}'.format(i)
                    if not os.path.exists(path):
                        os.mkdir(path)
                    tmp_path = path + '/tmp'
                    if not os.path.exists(tmp_path):
                        os.mkdir(tmp_path)
                outfq = args.output + '/{0}/{1}.fastq'.format(i, i)
                filter(args.input, outfq, i)
                sort_bam = args.output + '/{0}/{1}_sort.bam'.format(i, i)
                cmd = 'minimap2 -ax map-ont {ref} {fastq} -t {t1}|samtools sort -@ {t2} -O bam -o {bam}'.format(
                    ref=args.genome,
                    fastq=outfq,
                    t1=args.t,
                    t2=args.t,
                    bam=sort_bam)
                os.system(cmd)
                if args.q:
                    temp_bam = args.output + '/{0}/{1}_sort_mappingQ{q}.bam'.format(i, i, q=args.q)
                    cmd = 'samtools view {bam} -q {q} -b >{output_bam}'.format(bam=sort_bam, q=args.q,output_bam=temp_bam)
                    os.system(cmd)
                    sort_bam = temp_bam
                # 用pileup方法检测变异
                out_vcf = args.output + '/{0}/mpileup.vcf'.format(i)
                cmd = 'bcftools mpileup {bam} --fasta-ref {fa} --threads {t}> {vcf}'.format(bam=sort_bam, fa=args.genome,vcf=out_vcf,t = args.t)
                os.system(cmd)
                # bcftools call
                vcf = args.output + '/{0}/call.vcf'.format(i)
                cmd = 'bcftools call {vcf} -c -v -o {snp_vcf} --threads {t}'.format(vcf=out_vcf, snp_vcf=vcf,t=args.t)
                os.system(cmd)
                # filter
                if args.qual or args.minDP:
                    if args.qual and args.minDP:
                        file = 'Q{0}_DP{1}.vcf'.format(args.qual, args.minDP)
                        filter_vcf = args.output + '/{0}/'.format(i) + file
                    elif args.qual:
                        file = 'Q{0}.vcf'.format(args.qual)
                        filter_vcf = args.output + '/{0}/'.format(i) + file
                    elif args.minDP:
                        file = 'DP{0}.vcf'.format(args.minDP)
                        filter_vcf = args.output + '/{0}/'.format(i) + file
                else:
                    file = 'filter.vcf'
                    filter_vcf = args.output + '/{0}/'.format(i) + file
                cmd = 'bcftools filter -e \'QUAL<{0} | DP <{1}\' {2} > {3}'.format(args.qual, args.minDP, vcf,filter_vcf)
                os.system(cmd)
                # view
                snp_vcf = args.output + '/{0}/snp_'.format(i) + file
                snp_file = 'snp_'.format(i) + file
                cmd = 'bcftools view --include "TYPE==\'{0}\'" {1} > {2}'.format(args.type, filter_vcf, snp_vcf)
                os.system(cmd)
                # FQ0
                fq_vcf = args.output + '/{0}/snp_FQ0_'.format(i) + file
                cmd = 'bcftools filter -e \'FQ>0\' {0} > {1}'.format(snp_vcf, fq_vcf)
                os.system(cmd)
                # 注释
                id_vcf = args.output + '/{0}/snp_id_'.format(i) + file
                id_file = 'snp_id_'.format(i) + file
                cmd = 'java -jar {0} annotate {1} {2} >{3}'.format(SnpSift, vcf_annotate, fq_vcf, id_vcf)
                os.system(cmd)
                snpEff_vcf = args.output + '/{0}/snpEff_'.format(i) + file
                snpEff_file = 'snpEff_'.format(i) + file
                cmd = 'java -jar {0} -c {1} {2} {3} -v >{4}'.format(snpEff, snpEff_config,self.hg ,id_vcf, snpEff_vcf)
                os.system(cmd)
                filelist = os.listdir(path)
                for f in filelist:
                    try:
                        int(f)
                        filt_dir = True
                    except:
                        filt_dir = False
                    if f != snp_file and f and id_file and f != snpEff_file and f != 'tmp' and not filt_dir:
                        mv_file = args.output + '/{0}/'.format(i) + f
                        shutil.move(mv_file,tmp_path)
        else:
            tmp_path = args.output + '/tmp'
            if not os.path.exists(tmp_path):
                os.mkdir(tmp_path)
            sort_bam = args.output + '/sort.bam'
            cmd = 'minimap2 -ax map-ont {ref} {fastq} -t {t1}|samtools sort -@ {t2} -O bam -o {bam}'.format(
                ref=args.genome, fastq=args.input, t1=args.t, t2=args.t, bam=sort_bam)
            os.system(cmd)
            if args.q:
                temp_bam = args.output + '/sort_mappingQ{q}.bam'.format(q=args.q)
                cmd = 'samtools view {bam} -q {q} -b >{output_bam}'.format(bam=sort_bam, q=args.q, output_bam=temp_bam)
                os.system(cmd)
                sort_bam = temp_bam
            out_vcf = args.output + '/mpileup.vcf'
            cmd = 'bcftools mpileup {bam} --fasta-ref {fa} --threads {t}> {vcf}'.format(bam=sort_bam, fa=args.genome, vcf=out_vcf , t=args.t)
            os.system(cmd)
            # bcftools call
            vcf = args.output + '/call.vcf'
            cmd = 'bcftools call {vcf} -c -v -o {snp_vcf} --threads {t}'.format(vcf=out_vcf, snp_vcf=vcf ,t = args.t)
            os.system(cmd)
            # filter
            if args.qual or args.minDP:
                if args.qual and args.minDP:
                    file = 'Q{0}_DP{1}.vcf'.format(args.qual, args.minDP)
                    filter_vcf = args.output + '/' + file
                elif args.qual:
                    file = 'Q{0}.vcf'.format(args.qual)
                    filter_vcf = args.output + '/' + file
                elif args.minDP:
                    file = 'DP{0}.vcf'.format(args.minDP)
                    filter_vcf = args.output + '/' + file
            else:
                file = 'filter.vcf'
                filter_vcf = args.output + '/' + file
            cmd = 'bcftools filter -e \'QUAL<{0} | DP <{1}\' {2} > {3}'.format(args.qual, args.minDP, vcf, filter_vcf)
            os.system(cmd)
            # view
            snp_vcf = args.output + '/snp_' + file
            snp_file = 'snp_' + file
            cmd = 'bcftools view --include "TYPE==\'{0}\'" {1} > {2}'.format(args.type, filter_vcf, snp_vcf)
            os.system(cmd)
            # FQ0
            fq_vcf = args.output + '/snp_FQ0_' + file
            cmd = 'bcftools filter -e \'FQ>0\' {0} > {1}'.format(snp_vcf, fq_vcf)
            os.system(cmd)
            # 注释
            id_vcf = args.output + '/snp_id_' + file
            id_file = 'snp_id_' + file
            cmd = 'java -jar {0} annotate {1} {2} >{3}'.format(SnpSift, vcf_annotate, fq_vcf, id_vcf)
            os.system(cmd)
            snpEff_vcf = args.output + '/snpEff_' + file
            snpEff_file = 'snpEff_' + file
            cmd = 'java -jar {0} -c {1} {2} {3} -v >{4}'.format(snpEff, snpEff_config,self.hg ,id_vcf, snpEff_vcf)
            os.system(cmd)
            filelist = os.listdir(args.output)
            for f in filelist:
                try:
                    int(f)
                    filt_dir = True
                except:
                    filt_dir = False
                if f != snp_file and f and id_file and f != snpEff_file and f != 'tmp' and not filt_dir:
                    mv_file = args.output + '/' + f
                    shutil.move(mv_file, tmp_path)
class nanosv(object):
    def __init__(self ,parser,argv, config):
        self.argv = argv
        self.args = ''
        self.parser = parser
        self.config = config
    def extract_param(self):
        parser = self.parser
        parser.add_argument('-i', '--input', required=True, help='input reads in fastq format')
        parser.add_argument('-o', '--output', required=True,
                            help="Name your project output directory.")
        parser.add_argument('-g', '--genome', required=True,
                            help="Reference genome sequence in FASTA format.")
        parser.add_argument('-t', type=int, default=16,
                            help='Number of parallel processes. default:16')
        parser.add_argument('-p', action='store_true', default=False,
                            help="Whether to use Porechop to  remove the adapter sequence."
                                 " default： False")
        parser.add_argument('-q', type=int, default=None,
                            help='minimum mapping quality of sub-reads. default:None')
        parser.add_argument('-l', '--minLength', type=int,
                            required=False, default=None, nargs='*',
                            help="Filter minimum length,Multiple values can be entered ,No filtering by default")

        parser.add_argument('-c','--config',default='',
                            help='The parameters of NanoSV. '
                                 'Give the full path to your own ini file.')
        parser.add_argument('-s','--sambamba',default='',
                            help='The parameters of NanoSV.'
                                 'Give the full path to the sambamba or samtools executable.')
        parser.add_argument('-b','--bed',default='',
                            help='The parameters of NanoSV.'
                                 'Give the full path to your own bed file, used for coverage depth calculations.')
        parser.add_argument('-f','--snp_file',help='The parameters of NanoSV.'
                                                   'Give full path to the SNP variant file for phasing. Supporting file formats: BED and VCF')


        self.args =  parser.parse_args(self.argv)
    def config_praram(self):
        args = self.args
        if not os.path.exists(args.output):
            os.mkdir(args.output)
    def log(self):
        args = self.args
        logging.info('\n\
        NanoSV analysis \n\
        input :                 {0}\n\
        output :                {1}\n\
        reference :             {2}\n\
        threads:                {3}\n\
        PoreChop:               {4}\n\
        minimum length :        {5}\n\
        Minimum Mapping Qual :  {6}'.format(args.input,
                                          args.output,
                                          args.genome,
                                          args.t,
                                          args.p,
                                          args.minLength,
                                          args.q))
    def run_cmd(self):
        args = self.args
        if args.output[-1] == '/':
            args.output = args.output[:-1]
        if args.p:
            output_porechop = args.output + '/trimmed_adapter.fastq'
            cmd = 'porechop -t {t} -i {input} -o {output}'.format(t=args.t,
                                                                  input=args.input,
                                                                  output=output_porechop)
            os.system(cmd)
            args.input = output_porechop
        if args.minLength:
            l = args.minLength
            for i in l:
                path = args.output + '/{0}'.format(i)
                if not os.path.exists(path):
                    os.mkdir(path)
                outfq = args.output + '/{0}/{1}.fastq'.format(i, i)
                filter(args.input, outfq, i)
                sort_bam = args.output + '/{0}/{1}_sort.bam'.format(i,i)
                cmd = 'minimap2 -ax map-ont {ref} {fastq} -t {t1}|samtools sort -@ {t2} -O bam -o {bam}'.format(
                    ref=args.genome,
                    fastq=args.input,
                    t1=args.t,
                    t2=args.t,
                    bam=sort_bam)
                os.system(cmd)
                if args.q:
                    temp_bam = args.output + '/{0}/{1}Q{2}.bam'.format(i,i,args.q)
                    cmd = 'samtools view {bam} -q {q} -b >{output_bam}'.format(bam=sort_bam, q=args.q,
                                                                               output_bam=temp_bam)
                    os.system(cmd)
                    sort_bam = temp_bam
                cmd = 'samtools index {0}'.format(sort_bam)
                os.system(cmd)
                output_bw = args.output + '/{0}/{1}Q{2}coverage.bw'.format(i,i,args.q)
                cmd = 'bamCoverage -b {0} -o {1}'.format(sort_bam, output_bw)
                os.system(cmd)

                output_vcf = args.output + '/{0}/{1}Q{2}nanosv.vcf'.format(i,i,args.q)
                # cmd = 'NanoSV -t {0} -o {1} {2} -c {3}'.format(args.t, output_vcf, sort_bam, config_path)
                cmd = f'NanoSV -t {args.t} -o {output_vcf} {sort_bam}'
                if args.sambamba:
                    cmd = f'{cmd} -s {args.sambamba}'
                if args.config:
                    cmd = f'{cmd} -c {args.config}'
                if args.bed:
                    cmd = f'{cmd} -b {args.bed}'
                if args.snp_file:
                    cmd = f'{cmd} -f {args.snp_file}'
                os.system(cmd)
        else:
            sort_bam = args.output + '/sort_mapping.bam'
            cmd = 'minimap2 -ax map-ont {ref} {fastq} -t {t1}|samtools sort -@ {t2} -O bam -o {bam}'.format(ref=args.genome,
                                                                                                            fastq=args.input,
                                                                                                            t1=args.t,
                                                                                                            t2=args.t,
                                                                                                            bam=sort_bam)
            os.system(cmd)
            if args.q:
                temp_bam = args.output + '/sort_mappingQ{q}.bam'.format(q=args.q)
                cmd = 'samtools view {bam} -q {q} -b >{output_bam}'.format(bam=sort_bam, q=args.q, output_bam=temp_bam)
                os.system(cmd)
                sort_bam = temp_bam
            cmd = 'samtools index {0}'.format(sort_bam)
            os.system(cmd)
            output_bw = args.output + '/coverage.bw'
            cmd = 'bamCoverage -b {0} -o {1}'.format(sort_bam, output_bw)
            os.system(cmd)
            output_vcf = args.output + '/nanosv.vcf'
            # cmd = 'NanoSV -t {0} -o {1} {2} -c {3}'.format(args.t, output_vcf, sort_bam, config_path)
            cmd = f'NanoSV -t {args.t} -o {output_vcf} {sort_bam}'
            if args.sambamba:
                cmd = f'{cmd} -s {args.sambamba}'
            if args.config:
                cmd = f'{cmd} -c {args.config}'
            else:
                cmd = f'{cmd} -c {config_path}'
            if args.bed:
                cmd = f'{cmd} -b {args.bed}'
            if args.snp_file:
                cmd = f'{cmd} -f {args.snp_file}'
            os.system(cmd)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='cas_nanopore',
                                     usage='''cas_nanopore <subprogram> [options]
Commands:
megalodon               megalodon analysis
bcftools                bcftools analysis
nanosv                  nanosv analysis
''')
    subparsers = parser.add_subparsers()
    parser_megalodon = subparsers.add_parser(name = 'megalodon',description='Megalodon analysis',
                                            prog = 'cas_nanopore megalodon',
                                            usage = 'cas_nanopore megalodon [options]')
    parser_bcftools = subparsers.add_parser(name = 'bcftools',description='bcftools analysis',
                                            prog = 'cas_nanopore bcftools',
                                            usage = 'cas_nanopore bcftools [options]')
    parser_nanosv = subparsers.add_parser(name='NanoSV',description='NanoSV analysis',
                                            prog = 'cas_nanopore nanosv',
                                            usage = 'cas_nanopore nanosv [options]')
    logging.getLogger().setLevel(logging.INFO)
    config = configparser.ConfigParser()
    path1 =  f'{os.path.dirname(os.path.realpath(sys.executable))}/config.ini'
    path2 = f'{sys.path[0]}/config.ini'
    if os.path.exists(path1):
        config_path = path1
    elif os.path.exists(path2):
        config_path = path2
    config.read(config_path)
    if len(sys.argv) <= 1 :
        #必须使用三种方法一种
        parser.print_help()
        time.sleep(0.01)
        sys.stderr.write('\nNo argument given to cas9_nanopore'
                         '\nExiting\n')
        sys.exit(0)
    else:
        if sys.argv[1] == 'megalodon':
            procedure = megalodon(parser_megalodon,sys.argv[2:], config)
            procedure.extract_param()
            procedure.config_praram()
            procedure.run_cmd()
            procedure.run_logging()
        elif sys.argv[1] == 'bcftools':
            procedure = bcftools(parser_bcftools,sys.argv[2:], config)
            procedure.extract_param()
            procedure.config_praram()
            procedure.log()
            procedure.run_cmd()
        elif sys.argv[1] == 'nanosv':
            procedure = nanosv(parser_nanosv,sys.argv[2:], config)
            procedure.extract_param()
            procedure.config_praram()
            procedure.log()
            procedure.run_cmd()
        else:
            parser.print_help()