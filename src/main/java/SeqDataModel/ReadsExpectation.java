package SeqDataModel;

/**
 * Reads count expectation of IP or INPUT data for each genes in a individual.
 * for peak i in individual j, the reads expectation of IP and INPUT data calculate with following equations,
 *
 *          expect_ip,i,j = size-factor_ip,j * exp_i * k * p_i
 *
 *          expect_input,i,j = size-factor_input,j * exp_i
 *
 *  where size-factor_ip,j and size-factor_input,j are the j th sample IP and INPUT data size factor
 *        exp_i is the background expression in peak region
 *        k is the expansion effect of antigen
 *        p_i denotes as the methylation level
 */
public class ReadsExpectation {
    private int individualNumber, geneNumber;
    private int[][] ipGenesReads, inputGeneReads;
    private double expansionEffect;
    private double[] methylationLevel, individualIPSizeFactors, individualINPUTSizeFactors, geneBkgExp = null;
    private BackgroundExpression bkgExp;

    /**
     * Constructor
     * @param ipGenesReads IP gene reads, shape individual number × gene number
     * @param inputGeneReads INPUT gene reads, shape individual number × gene number
     * @param methylationLevel methylation level, shape 1 × gene number
     * @param expansionEffect expansion effect of antigen
     */
    public ReadsExpectation(int[][] ipGenesReads, int[][] inputGeneReads,
                            double[] methylationLevel, double expansionEffect) {
        assert ipGenesReads.length == inputGeneReads.length;
        this.individualNumber = ipGenesReads.length;
        this.geneNumber = ipGenesReads[0].length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
        this.expansionEffect = expansionEffect;
        this.bkgExp = new BackgroundExpression(ipGenesReads, inputGeneReads);
        this.methylationLevel = methylationLevel;
    }

    /**
     * Constructor
     * @param ipGenesReads IP gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param inputGeneReads INPUT gene reads, shape m×n, where m denotes as individual number, n denote as gene number
     * @param methylationLevel methylation level, shape 1×n, where n represents gene number
     */
    public ReadsExpectation(int[][] ipGenesReads, int[][] inputGeneReads,
                            double[] methylationLevel, double[] bkgExp) {
        assert ipGenesReads.length == inputGeneReads.length;
        this.individualNumber = ipGenesReads.length;
        this.geneNumber = ipGenesReads[0].length;
        this.ipGenesReads = ipGenesReads;
        this.inputGeneReads = inputGeneReads;
        this.bkgExp = new BackgroundExpression(ipGenesReads, inputGeneReads);
        this.methylationLevel = methylationLevel;
        this.geneBkgExp = bkgExp;
    }

    /**
     * use median of total genes' size factor as size factor for each IP or INPUT individual
     * @param input true, if get INPUT data size factor; otherwise, IP data size factor
     */
    private void globalSizeFactor(boolean input) {
        // shape 1 × individualNumber
        double[] individualSizeFactors = new double[this.individualNumber];

        for (int i=0; i<this.individualNumber; i++) {
            // reads count of genes in IP and INPUT data of individual i, shape 1 × geneNumber
            int[] individualIPReads = this.ipGenesReads[i], individualINPUTReads = this.inputGeneReads[i];
            SizeFactor sf = new SizeFactor(individualIPReads, individualINPUTReads);
            double sizeFactor = sf.getSizeFactor(input);
            sf = null;
            individualSizeFactors[i] = sizeFactor;
            // only for simulation data test
            // individualSizeFactors[i] = 1;
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
     * calculate gene expression for each gene with specific size factor
     */
    private void calcGenesBkgExp(double[] sizeFactors) {
        this.geneBkgExp = this.bkgExp.geneBackgroundExp(sizeFactors);
    }

    /**
     * calculate peak region reads count expectations for each gene in individuals
     * @param input true, if get INPUT data reads count expectation; otherwise, IP data reads count expectation
     * @return reads count expectation array. In shape m×n, where m denotes as individual number, n denote as gene number
     */
    private double[][] geneReadsExpectation(boolean input) {
        double[] individualSizeFactors = input? this.individualINPUTSizeFactors: this.individualIPSizeFactors;

        double[][] readsExpectation = new double[this.individualNumber][this.geneNumber];

        for (int i=0; i<this.individualNumber; i++) {
            double individualSizeFactor = individualSizeFactors[i];

            for (int j=0; j< this.geneNumber; j++) {
                if (!input)
                    individualSizeFactor = individualSizeFactor * this.methylationLevel[j] * this.expansionEffect;

                double bkgExp = this.geneBkgExp[j];  // background expression of the i th peak region
                readsExpectation[i][j] = individualSizeFactor * bkgExp;
            }
        }
        return readsExpectation;
    }

    /**
     * genes' reads count expectation in IP data of each individual
     * @return genes' reads count expectation
     */
    public double[][] getIPReadsExpectation() {
        if (this.geneBkgExp == null)
            this.calcGenesBkgExp();
        this.globalSizeFactor(false);
        // only for simulation data test
        // this.individualIPSizeFactors = new double[] {1, 1, 1};
        return this.geneReadsExpectation(false);
    }

    /**
     * genes' reads count expectation in IP data of each individual with specific size factors
     * @return genes' reads count expectation
     */
    public double[][] getIPReadsExpectation(double[] individualIPSizeFactors) {
        this.calcGenesBkgExp(individualIPSizeFactors);
        this.individualIPSizeFactors = individualIPSizeFactors;
        // only for simulation data test
        // this.individualIPSizeFactors = new double[] {1, 1, 1};
        return this.geneReadsExpectation(false);
    }

    /**
     * genes' reads count expectation in INPUT data of each individual
     * @return genes' reads count expectation
     */
    public double[][] getINPUTReadsExpectation() {
        if (this.geneBkgExp == null)
            this.calcGenesBkgExp();
        this.globalSizeFactor(true);
        // only for simulation data test
        // this.individualINPUTSizeFactors = new double[] {1, 1, 1};
        return this.geneReadsExpectation(true);
    }

    /**
     * genes' reads count expectation in INPUT data of each individual with specific size factors
     * @return genes' reads count expectation
     */
    public double[][] getINPUTReadsExpectation(double[] individualINPUTSizeFactors) {
        this.calcGenesBkgExp(individualINPUTSizeFactors);
        this.individualINPUTSizeFactors = individualINPUTSizeFactors;
        // only for simulation data test
        // this.individualINPUTSizeFactors = new double[] {1, 1, 1};
        return this.geneReadsExpectation(true);
    }
}
