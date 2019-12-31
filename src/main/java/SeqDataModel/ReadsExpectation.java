package SeqDataModel;

/**
 * Reads count expectation of IP or INPUT data for each genes in a individual
 */
public class ReadsExpectation {
    private int[][] ipGenesReads, inputGeneReads;
    private double[] methylationLevel, individualIPSizeFactors, individualINPUTSizeFactors, geneBkgExp = null;
    private BackgroundExpression bkgExp;

    /**
     * Constructor
     * @param ipGenesReads IP gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param inputGeneReads INPUT gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param methylationLevel methylation level, shape 1×n, where n represents gene number
     */
    public ReadsExpectation(int[][] ipGenesReads, int[][] inputGeneReads, double[] methylationLevel) {
        assert ipGenesReads.length == inputGeneReads.length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
        this.bkgExp = new BackgroundExpression(ipGenesReads, inputGeneReads);
        this.methylationLevel = methylationLevel;
    }

    /**
     * Constructor
     * @param ipGenesReads IP gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param inputGeneReads INPUT gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param methylationLevel methylation level, shape 1×n, where n represents gene number
     * @param geneBkgExp gene background expression, shape 1×n, where n represents gene number
     */
    public ReadsExpectation(int[][] ipGenesReads, int[][] inputGeneReads, double[] methylationLevel, double[] geneBkgExp) {
        assert ipGenesReads.length == inputGeneReads.length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
        this.bkgExp = new BackgroundExpression(ipGenesReads, inputGeneReads);
        this.methylationLevel = methylationLevel;
        this.geneBkgExp = geneBkgExp;
    }

    /**
     * use median of total genes' size factor as size factor for each IP or INPUT individual
     * @param input true, if get INPUT data size factor; otherwise, IP data size factor
     */
    private void globalSizeFactor(boolean input) {
        // shape 1 × individualNumber
        double[] individualSizeFactors = new double[ipGenesReads.length];

        for (int i=0; i<individualSizeFactors.length; i++) {
            // reads count of genes in IP and INPUT data of individual i, shape 1 × geneNumber
            int[] individualIPReads = this.ipGenesReads[i], individualINPUTReads = this.inputGeneReads[i];
            SizeFactor sf = new SizeFactor(individualIPReads, individualINPUTReads);
            double sizeFactor = sf.getSizeFactor(input);
            sf = null;
            individualSizeFactors[i] = sizeFactor;
        }

        if (input)
            this.individualINPUTSizeFactors = individualSizeFactors;
        else
            this.individualIPSizeFactors = individualSizeFactors;
    }

    /**
     * calculate gene expression for each gene
     */
    private void calcGenesBkgExp() {
        this.geneBkgExp = this.bkgExp.geneBackgroundExp();
    }

    /**
     * calculate reads count expectations for each gene in individuals
     * @param input true, if get INPUT data reads count expectation; otherwise, IP data reads count expectation
     * @return reads count expectation array. In shape m×n, where m denotes as individual number, n denote as gene number
     */
    private double[][] geneReadsExpectation(boolean input) {
        double[] individualSizeFactors = input? this.individualINPUTSizeFactors: this.individualIPSizeFactors;

        int individualNum = individualSizeFactors.length, geneNum = this.geneBkgExp.length;
        double[][] readsExpectation = new double[individualNum][geneNum];

        for (int i=0; i<individualNum; i++) {
            double individualSizeFactor = individualSizeFactors[i];

            for (int j=0; j< geneNum; j++) {
                if (!input)
                    individualSizeFactor = individualSizeFactor * this.methylationLevel[j];
                double bkgExp = this.geneBkgExp[j];
                readsExpectation[i][j] = individualSizeFactor * bkgExp;
            }
        }
        return readsExpectation;
    }

    /**
     * genes' reads count expectation in IP data of each individual
     * @return genes' reads count expectation
     */
    public double[][] getIPReadsExepectation() {
        if (this.geneBkgExp == null)
            this.calcGenesBkgExp();
        this.globalSizeFactor(false);
//        // TODO:当作用于实际数据时改用上一行被注释的代码
//        this.individualIPSizeFactors = new double[] {1, 1, 1};
        return this.geneReadsExpectation(false);
    }

    public double[][] getINPUTReadsExpectation() {
        if (this.geneBkgExp == null)
            this.calcGenesBkgExp();
        this.globalSizeFactor(true);
        // TODO:当作用于实际数据时改用上一行被注释的代码
//        this.individualINPUTSizeFactors = new double[] {1, 1, 1};
        return this.geneReadsExpectation(true);
    }
}
