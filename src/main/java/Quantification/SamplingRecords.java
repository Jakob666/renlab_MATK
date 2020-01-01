package Quantification;

public class SamplingRecords {
    private int[][] sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads;
    private double geneBackgroundExp, nonPeakExpression;
    private double[] nonSpecificEnrichment, methylationLevels, ipOverdispersions, inputOverdispersions, backgroundExpressions;
    private double[][] sampleIPExpectations, sampleINPUTExpectations, nonPeakIPExpectations, nonPeakINPUTExpectations;

    public SamplingRecords() {}

    public void setReads(int[][] sampleIPReads, int[][] sampleINPUTReads,
                         int[][] nonPeakIPReads, int[][] nonPeakINPUTReads) {
        this.sampleIPReads = sampleIPReads;
        this.sampleINPUTReads = sampleINPUTReads;
        this.nonPeakIPReads = nonPeakIPReads;
        this.nonPeakINPUTReads = nonPeakINPUTReads;
    }

    public void setExpectations(double[][] sampleIPExpectations, double[][] sampleINPUTExpectations,
                                double[][] nonPeakIPExpectations, double[][] nonPeakINPUTExpectations) {
        this.sampleIPExpectations = sampleIPExpectations;
        this.sampleINPUTExpectations = sampleINPUTExpectations;
        this.nonPeakIPExpectations = nonPeakIPExpectations;
        this.nonPeakINPUTExpectations = nonPeakINPUTExpectations;
    }

    public void setSampleIPExpectations(double[][] sampleIPExpectations) {
        this.sampleIPExpectations = sampleIPExpectations;
    }

    public void setNonPeakIPExpectations(double[][] nonPeakIPExpectations) {
        this.nonPeakIPExpectations = nonPeakIPExpectations;
    }

    public int[][] getSampleIPReads() {
        return this.sampleIPReads;
    }

    public int[][] getSampleINPUTReads() {
        return this.sampleINPUTReads;
    }

    public int[][] getNonPeakIPReads() {
        return this.nonPeakIPReads;
    }

    public int[][] getNonPeakINPUTReads() {
        return this.nonPeakINPUTReads;
    }

    public double[][] getSampleIPExpectations() {
        return this.sampleIPExpectations;
    }

    public double[][] getSampleINPUTExpectations() {
        return this.sampleINPUTExpectations;
    }

    public double[][] getNonPeakIPExpectations() {
        return this.nonPeakIPExpectations;
    }

    public double[][] getNonPeakINPUTExpectations() {
        return this.nonPeakINPUTExpectations;
    }

    public void setNonSpecificEnrichment(double[] nonSpecificEnrichment) {
        this.nonSpecificEnrichment = nonSpecificEnrichment;
    }

    public double getNonSpecificEnrichmentRatio(int time) {
        return this.nonSpecificEnrichment[time];
    }

    public void setNonSpecificEnrichmentRatio(double value, int time) {
        this.nonSpecificEnrichment[time] = value;
    }

    public double[] getNonSpecificEnrichment() {
        return this.nonSpecificEnrichment;
    }

    public void setMethylationLevels(double[] methylationLevels) {
        this.methylationLevels = methylationLevels;
    }

    public void setMethylationLevelValue(double value, int time) {
        this.methylationLevels[time] = value;
    }

    public double getMethylationLevelValue(int time) {
        return this.methylationLevels[time];
    }

    public double[] getMethylationLevels() {
        return this.methylationLevels;
    }

    public void setIpOverdispersions(double[] ipOverdispersions) {
        this.ipOverdispersions = ipOverdispersions;
    }

    public void setIpOverdispersionValue(double value, int time) {
        this.ipOverdispersions[time] = value;
    }

    public double getIpOverdispersionValue(int time) {
        return this.ipOverdispersions[time];
    }

    public double[] getIpOverdispersions() {
        return this.ipOverdispersions;
    }

    public void setInputOverdispersions(double[] inputOverdispersions) {
        this.inputOverdispersions = inputOverdispersions;
    }

    public void setInputOverdispersionValue(double value, int time) {
        this.inputOverdispersions[time] = value;
    }

    public double getInputOverdispersionValue(int time) {
        return this.inputOverdispersions[time];
    }

    public double[] getInputOverdispersions() {
        return this.inputOverdispersions;
    }

    public void setBackgroundExpressions(double[] backgroundExpressions) {
        this.backgroundExpressions = backgroundExpressions;
    }

    public double getBackgroundExpressionValue(int time) {
        return this.backgroundExpressions[time];
    }

    public double[] getBackgroundExpressions() {
        return this.backgroundExpressions;
    }

    public void setGeneBackgroundExp(double geneBackgroundExp) {
        this.geneBackgroundExp = geneBackgroundExp;
    }

    public double getGeneBackgroundExp() {
        return this.geneBackgroundExp;
    }

    public void setNonPeakExpression(double nonPeakExpression) {
        this.nonPeakExpression = nonPeakExpression;
    }

    public double getNonPeakExpression() {
        return this.nonPeakExpression;
    }
}
