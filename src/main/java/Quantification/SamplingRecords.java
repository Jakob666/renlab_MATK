package Quantification;

public class SamplingRecords {
    private int[][] sampleIPReads, sampleINPUTReads, nonPeakIPReads, nonPeakINPUTReads;
    private double[] nonSpecificEnrichment, methylationLevels, ipOverdispersions, inputOverdispersions, backgroundExpressions, nonPeakExpressions;
    private double[] ipSizeFactors, inputSizeFactors, ipNonPeakSizeFactors, inputNonPeakSizeFactors;
    private double[][] sampleIPExpectations, sampleINPUTExpectations, nonPeakIPExpectations, nonPeakINPUTExpectations;
    private BackgroundExpressionSampler bkgExpSampler, nonPeakExpSampler;

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

    public void setSampleINPUTExpectations(double[][] sampleINPUTExpectations) {
        this.sampleINPUTExpectations = sampleINPUTExpectations;
    }

    public void setNonPeakIPExpectations(double[][] nonPeakIPExpectations) {
        this.nonPeakIPExpectations = nonPeakIPExpectations;
    }

    public void setNonPeakINPUTExpectations(double[][] nonPeakINPUTExpectations) {
        this.nonPeakINPUTExpectations = nonPeakINPUTExpectations;
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

    public void setBkgExpSampler(BackgroundExpressionSampler sampler) {
        this.bkgExpSampler = sampler;
    }

    public BackgroundExpressionSampler getBkgExpSampler() {
        return this.bkgExpSampler;
    }

    public void setBackgroundExpressions(double[] backgroundExpressions) {
        this.backgroundExpressions = backgroundExpressions;
    }

    public void setBackgroundExpressionValue(double value, int time) {
        this.backgroundExpressions[time] = value;
    }

    public double getBackgroundExpressionValue(int time) {
        return this.backgroundExpressions[time];
    }

    public double[] getBackgroundExpressions() {
        return this.backgroundExpressions;
    }

    public void setNonPeakExpSampler(BackgroundExpressionSampler sampler) {
        this.nonPeakExpSampler = sampler;
    }

    public BackgroundExpressionSampler getNonPeakExpSampler() {
        return this.nonPeakExpSampler;
    }

    public void setNonPeakExpressions(double[] nonPeakExpressions) {
        this.nonPeakExpressions = nonPeakExpressions;
    }

    public void setNonPeakExpressionValue(double value, int time) {
        this.nonPeakExpressions[time] = value;
    }

    public double getNonPeakExpressionValue(int time) {
        return this.nonPeakExpressions[time];
    }

    public double[] getNonPeakExpressions() {
        return this.nonPeakExpressions;
    }

    public void setInputSizeFactors(double[] inputSizeFactors) {
        this.inputSizeFactors = inputSizeFactors;
    }

    public double[] getInputSizeFactors() {
        return this.inputSizeFactors;
    }

    public void setInputNonPeakSizeFactors(double[] inputNonPeakSizeFactors) {
        this.inputNonPeakSizeFactors = inputNonPeakSizeFactors;
    }

    public double[] getInputNonPeakSizeFactors() {
        return this.inputNonPeakSizeFactors;
    }

    public void setIpSizeFactors(double[] ipSizeFactors) {
        this.ipSizeFactors = ipSizeFactors;
    }

    public double[] getIpSizeFactors() {
        return this.ipSizeFactors;
    }

    public void setIpNonPeakSizeFactors(double[] ipNonPeakSizeFactors) {
        this.ipNonPeakSizeFactors = ipNonPeakSizeFactors;
    }

    public double[] getIpNonPeakSizeFactors() {
        return this.ipNonPeakSizeFactors;
    }
}
