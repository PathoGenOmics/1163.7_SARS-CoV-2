#! /usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import subprocess
def parse_args():
    '''
    Parse arguments given to script
    '''
    parser = argparse.ArgumentParser(description = "Parse COVs to get mutations")
    parser.add_argument("-vcf_file", dest = "vcf_file", required = True)
    parser.parse_args()
    args = parser.parse_args()

    return args

def vcfmaker(args):
    vcf_file = args.vcf_file
    with open(vcf_file+".freq","w") as out_freqaln:
        with open(vcf_file,'r') as input_aln:
            tit_inicial="#Position\tTotalMut\tFrequency\tPercentage\tVariants\n"
            out_freqaln.write(tit_inicial)
            for line in input_aln:
                if "#" not in line:
                    linea=line.replace("\n","").split("\t")
                    listanum=linea[9:]
                    #print(listanum)
                    alternativa = linea[4].split(",")
                    total=0
                    sentencia=''
                    for i in range(len(alternativa)):
                        if alternativa[i] != "*":
                            numero = i+1
                            cont_alt = listanum.count(str(numero))
                            total += cont_alt
                            sentencia += alternativa[i]+":"+str(cont_alt)+", "
                    frase=linea[1]+"\t"+str(total)+"\t"+str(total/(len(listanum)-1))+"\t"+str(total/(len(listanum)-1)*100)+"\t"+sentencia+"\n"
                    out_freqaln.write(frase)    

def anotarSnpEff(args):
    vcf_file = args.vcf_file
    with open('outMut'+ vcf_file + '.parannot','w') as out_freqaln:
        with open(vcf_file,'r') as input_aln:
            for line in input_aln:
                if "##" in line:
                    out_freqaln.write(line)
                elif "#" in line:
                    line=line.replace("\n","").split("\t")
                    titulillo=line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\n"
                    out_freqaln.write(titulillo)
                else:
                    line=line.replace("\n","").split("\t")
                    linea="NC_045512.2\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]+"\t"+line[6]+"\t"+line[7]+"\n"
                    out_freqaln.write(linea)                    
        subprocess.call('java -Xmx4g -jar /home/snpEff/snpEff/snpEff.jar ann -c /home/snpEff/snpEff/snpEff.config -noStats -no-downstream -no-upstream NC_045512.2 '+ 'outMut'+ today + '.parannot >'+ today + 'res.annot', shell = True)
        subprocess.call('rm outMut'+ vcf_file + '.parannot ', shell = True)

def main():
    args = parse_args() #Recogemos los argumentos
    vcfmaker(args)
    anotarSnpEff(args)
    
if __name__ == '__main__':
    main()