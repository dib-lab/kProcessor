kProcessor::loadFromKMC(currentFrame,filename);
// calculate the sum of all kmers for normalization
any totalCountAny=kProcessor::aggregate(currentFrame,(uint64_t)0,  [](kmerRow it, any v) -> any {
    uint32_t count0;
    it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count",count0);
    return (any)(any_cast<uint64_t>(v) + (uint64_t)count0);
});
// normalize the counts kmers
currentFrame= kProcessor::transform(currentFrame,  [=](kmerRow it) -> kmerRow {
    uint32_t count0;
    it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count",count0);
    double normalized = (double)count0*(100000000.0) / totalCount;
    it.setColumnValue<uint32_t,vectorColumn<uint32_t> >("count",(uint32_t)normalized);
    return it;
});
// index reference
kProcessor::index(KMERS, genes_file+".names", genesFrame);


// join all kdataframes. required Indeices contains the index of the index frame only
// remove kmers that has zero counts in all samples
kDataFrame* res= kProcessor::innerJoin(kFrames, requiredIndices);
res=kProcessor::filter(res,[=](kmerRow r) -> bool {
    for(unsigned i=0; i < allDatasets ;i++ ){
        uint32_t count;
        r.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count.0",count);
        if(count>0)
            return true;
    }
    return false;
});


// calculate fold change
res->addColumn("foldChange",new vectorColumn<double >(res->size()));
res= kProcessor::transform(res,  [=](kmerRow it) -> kmerRow {
    uint32_t sample_count0, sample_count1, sample_count2;
    uint32_t control_count0, control_count1, control_count2;
    it.getColumnValue<uint32_t,vectorColumn<uint32_t> >("count.0",sample_count0);
    // get all counts similar to the above line

    double sampleAVG= (double)(sample_count0 + sample_count1 + sample_count2) / (double)nSamples;
    double controlAVG= (double)(control_count0 + control_count1 + control_count2) / (double)nControl;
    double foldChange= sampleAVG / controlAVG;
    it.setColumnValue<double,vectorColumn<double> >("foldChange",foldChange);
    return it;
});


// gather fold changes by genes
auto foldChangeByGene=new unordered_map<uint32_t ,vector<double> >();
any genesGatherAny=kProcessor::aggregate(res,foldChangeByGene,  [=](kmerRow it, any v) -> any {
    auto dict=any_cast<unordered_map<uint32_t ,vector<double>>*>(v);
    double foldChange;
    it.getColumnValue<double,vectorColumn<double> >("foldChange",foldChange);
    vector<uint32_t> color;
    it.getColumnValue<vector<uint32_t> ,
            deduplicatedColumn<vector<uint32_t>,StringColorColumn> >(colorColumn,color);
    for(auto c: color)
    {
        (*dict)[c].push_back(foldChange);
    }
    return (any)(dict);
});