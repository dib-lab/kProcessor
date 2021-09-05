import kProcessor as kp
import faulthandler
faulthandler.enable()

def differntialExpression(genes_file, samplesInput, controlInput, outputFilename):
    nSamples = len(samplesInput)
    nControl = len(controlInput)
    allDatasets = nSamples + nControl

    kFrames = list()  # kp.kFramesVector
    allSamples = list()
    allSamples = samplesInput + controlInput

    requiredIndices = list()

    index = 0

    count_col_name = "count"

    for filename in allSamples:
        currentFrame = kp.kDataFrameMAP()
        print(f"Loading {filename} KMC DB")
        kp.loadFromKMC(currentFrame, filename)
        print(f"Loading {filename} kmers: {currentFrame.size()}")
        totalCount = kp.aggregate_count(currentFrame, count_col_name)
        print(f"Total count = {totalCount}")
        kp.transform_normalize(currentFrame, count_col_name, totalCount)
        kFrames.append(currentFrame)

    kSize = int()
    if len(kFrames):
        kSize = kFrames[-1].getkSize()
        print(f"[DEBUG] kSize = {kSize}")

    chunkSize = 1000
    genesFrame = kp.kDataFrameMAP(kSize)
    kp.index(genesFrame, {"kSize": kSize}, genes_file, chunkSize, f"{genes_file}.names")
    # kp.createColorColumn(genesFrame)
    kFrames.append(genesFrame)
    requiredIndices.append(len(kFrames)-1)
    colorColumn = f"color{len(kFrames)-1}"
    print(f"Load {genes_file} kmers: {kFrames[-1].size()}")
    res = kp.innerJoin(kFrames, requiredIndices)
    print(f"Joined {res.size()} kmers...")
    print("[DEBUG] filtering zeroCounts...")
    res = kp.filter_zeroCounts(res, allDatasets)
    foldChange_col_name = "foldChange"
    doubleVectorColumn = kp.vectorColumn_double(res.size())
    res.addColumn(foldChange_col_name, doubleVectorColumn)
    print("[DEBUG] transforming foldChanges...")
    kp.transform_foldchange(res, nSamples, nControl, allDatasets, foldChange_col_name)
    print("[DEBUG] aggregating foldChangeByGene...")

    foldChangeByGene = kp.aggregate_foldChangeByGene(res, foldChange_col_name, colorColumn)

    it = foldChangeByGene.begin()
    with open(output_filename, 'w') as output:
        while(it != foldChangeByGene.end()):
            k, v = it.next()
            if len(v):
                v = sorted(v)
                median = v[len(v) // 2]
                output.write(f"{k}\t{median}\n")


# Inputs
input_genes_file = "genesERCC.fa" #str()
input_KMC_sample_files = ["sample.kmc.db"] # list()
input_KMC_control_files = ["control.kmc.db"] #list()
output_filename = "pykDiff_sample_control" #str()

differntialExpression(input_genes_file, input_KMC_sample_files, input_KMC_control_files, output_filename)
