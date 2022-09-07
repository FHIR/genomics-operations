import os
import subprocess
import sys

fileFolderPath = "C:\\snpEff\\simple\\38"
outputFolderPath = "C:\\snpEff\\output"
snpEffPath = "C:\\snpEff"
databasePath = "C:\\snpEff\\gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz"

genome = 'GRCh38.p13'

snpEffJarPath = os.path.join(snpEffPath, 'snpEff.jar')
snpSiftJarPath = os.path.join(snpEffPath, 'snpSift.jar')


if sys.platform == "Windows":
    MOVE = "move"
    MKDIR = "mkdir"
else:
    MOVE = "mv"
    MKDIR = "md"

os.chdir(snpEffPath)

for vcfFile in os.listdir(fileFolderPath):
    print('in file: ' + vcfFile)
    filePath = fileFolderPath + '\\' + vcfFile
    outputPath = outputFolderPath + '\\' + vcfFile
    snpEffCmd = ['java', '-jar', snpEffJarPath] + ['-o', 'vcf'] + [genome, filePath]

    intVCFPath = os.path.join(snpEffPath, 'intVCF.vcf')
    with open(intVCFPath, 'w') as intVCF:
        subprocess.run(snpEffCmd, stdout=intVCF)
        intVCF.close()

    changedAFPath = os.path.join(snpEffPath, 'changedAF.vcf')
    with open(intVCFPath, "r") as orig:
        with open(changedAFPath, 'w') as updated:
            for line in orig:
                newString = line.replace(';AF=', ';CHANGEDAF=')
                updated.write(newString)
            
    snpSiftCmd = 'java -jar ' + snpSiftJarPath + ' annotate ' + databasePath + ' -info "AF" ' + changedAFPath

    annotatedVCFPath = os.path.join(snpEffPath, 'annotatedVCF.vcf')

    with open(annotatedVCFPath, 'w') as annotatedFile:
        subprocess.run(snpSiftCmd, stdout=annotatedFile)

    with open(annotatedVCFPath, "rt") as orig:
        with open(outputPath, "wt") as updated:
            for line in orig:
                newLine = line.replace(';AF=', ';POPAF=')
                newLine = newLine.replace(';CHANGEDAF=', ';AF=')
                updated.write(newLine)
            
    os.remove(intVCFPath)
    os.remove(annotatedVCFPath)
    os.remove(changedAFPath)
    