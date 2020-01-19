package Quantification;

public class SamplingRecords {
    private double[] methylationLevels, ipOverdispersions, inputOverdispersions, backgroundExpressions, expansionEffects;
    private double[] ipSizeFactors, inputSizeFactors;

    public SamplingRecords() {}

    public void setMethylationLevels(double[] methylationLevels) {
        this.methylationLevels = methylationLevels;
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

    public void setExpansionEffects(double[] expansionEffects) {
        this.expansionEffects = expansionEffects;
    }

    public void setExpansionEffectsValue(double value, int time) {
        this.expansionEffects[time] = value;
    }

    public double getExpansionEffectsValue(int time) {
        return this.expansionEffects[time];
    }

    public double[] getExpansionEffects() {
        return this.expansionEffects;
    }

    public void setBackgroundExpressions(double[] backgroundExpressions) {
        this.backgroundExpressions = backgroundExpressions;
    }

    public double[] getBackgroundExpressions() {
        return this.backgroundExpressions;
    }

    public void setInputSizeFactors(double[] inputSizeFactors) {
        this.inputSizeFactors = inputSizeFactors;
    }

    public double[] getInputSizeFactors() {
        return this.inputSizeFactors;
    }

    public void setIpSizeFactors(double[] ipSizeFactors) {
        this.ipSizeFactors = ipSizeFactors;
    }

    public double[] getIpSizeFactors() {
        return this.ipSizeFactors;
    }
}
