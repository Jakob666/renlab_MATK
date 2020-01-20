package DifferentialMethylation;

public class Parameters {
    private double quantifiedTretIPOverdispersion, quantifiedTretINPUTOverdispersion
                   , quantifiedCtrlIPOverdispersion, quantifiedCtrlINPUTOverdispersion, quantifiedExpansionEffect;
    private double[] quantifiedTretMethLevel, quantifiedCtrlMethLevel, quantifiedTretExpression, quantifiedCtrlExpression;
    // shape 1 × samplingTime
    private double[] tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion;
    // shape 1 × peakNumber
    private double[] tretMethylation, ctrlMethylation, tretBkgExp, ctrlBkgExp, expansionEffect;
    // shape 1 × individualNumber
    private double[] tretIPSizeFactors, tretINPUTSizeFactors, ctrlIPSizeFactors, ctrlINPUTSizeFactors;

    public Parameters() {}

    public double[] getTretIPOverdispersion() {
        return this.tretIPOverdispersion;
    }

    public void setTretIPOverdispersion(double[] value) {
        this.tretIPOverdispersion = value;
    }

    public double[] getTretINPUTOverdispersion() {
        return this.tretINPUTOverdispersion;
    }

    public void setTretINPUTOverdispersion(double[] value) {
        this.tretINPUTOverdispersion = value;
    }

    public double[] getCtrlIPOverdispersion() {
        return this.ctrlIPOverdispersion;
    }

    public void setCtrlIPOverdispersion(double[] value) {
        this.ctrlIPOverdispersion = value;
    }

    public double[] getCtrlINPUTOverdispersion() {
        return this.ctrlINPUTOverdispersion;
    }

    public void setCtrlINPUTOverdispersion(double[] value) {
        this.ctrlINPUTOverdispersion = value;
    }

    public double[] getTretMethylation() {
        return this.tretMethylation;
    }

    public void setTretMethylation(double[] value) {
        this.tretMethylation = value;
    }

    public double[] getCtrlMethylation() {
        return this.ctrlMethylation;
    }

    public void setCtrlMethylation(double[] value) {
        this.ctrlMethylation = value;
    }

    public void setExpansionEffect(double[] value) {
        this.expansionEffect = value;
    }

    public double[] getExpansionEffect() {
        return this.expansionEffect;
    }

    public double[] getTretBkgExp() {
        return this.tretBkgExp;
    }

    public void setTretBkgExp(double[] value) {
        this.tretBkgExp = value;
    }

    public double[] getCtrlBkgExp() {
        return this.ctrlBkgExp;
    }

    public void setCtrlBkgExp(double[] value) {
        this.ctrlBkgExp = value;
    }

    public void setTretINPUTSizeFactors(double[] tretINPUTSizeFactors) {
        this.tretINPUTSizeFactors = tretINPUTSizeFactors;
    }

    public double[] getTretINPUTSizeFactors() {
        return this.tretINPUTSizeFactors;
    }

    public void setTretIPSizeFactors(double[] tretIPSizeFactors) {
        this.tretIPSizeFactors = tretIPSizeFactors;
    }

    public double[] getTretIPSizeFactors() {
        return this.tretIPSizeFactors;
    }

    public void setCtrlINPUTSizeFactors(double[] ctrlINPUTSizeFactors) {
        this.ctrlINPUTSizeFactors = ctrlINPUTSizeFactors;
    }

    public double[] getCtrlINPUTSizeFactors() {
        return this.ctrlINPUTSizeFactors;
    }

    public void setCtrlIPSizeFactors(double[] ctrlIPSizeFactors) {
        this.ctrlIPSizeFactors = ctrlIPSizeFactors;
    }

    public double[] getCtrlIPSizeFactors() {
        return this.ctrlIPSizeFactors;
    }

    public void setQuantifiedCtrlExpression(double[] quantifiedCtrlExpression) {
        this.quantifiedCtrlExpression = quantifiedCtrlExpression;
    }

    public double[] getQuantifiedCtrlExpression() {
        return this.quantifiedCtrlExpression;
    }

    public void setQuantifiedTretExpression(double[] quantifiedTretExpression) {
        this.quantifiedTretExpression = quantifiedTretExpression;
    }

    public double[] getQuantifiedTretExpression() {
        return this.quantifiedTretExpression;
    }

    public void setQuantifiedTretMethLevel(double[] quantifiedTretMethLevel) {
        this.quantifiedTretMethLevel = quantifiedTretMethLevel;
    }

    public double[] getQuantifiedTretMethLevel() {
        return this.quantifiedTretMethLevel;
    }

    public void setQuantifiedCtrlMethLevel(double[] quantifiedCtrlMethLevel) {
        this.quantifiedCtrlMethLevel = quantifiedCtrlMethLevel;
    }

    public double[] getQuantifiedCtrlMethLevel() {
        return this.quantifiedCtrlMethLevel;
    }

    public void setOverdispersions(double quantifiedTretIPOverdispersion, double quantifiedTretINPUTOverdispersion,
                                   double quantifiedCtrlIPOverdispersion, double quantifiedCtrlINPUTOverdispersion) {
        this.quantifiedTretIPOverdispersion = quantifiedTretIPOverdispersion;
        this.quantifiedCtrlIPOverdispersion = quantifiedCtrlIPOverdispersion;
        this.quantifiedTretINPUTOverdispersion = quantifiedTretINPUTOverdispersion;
        this.quantifiedCtrlINPUTOverdispersion = quantifiedCtrlINPUTOverdispersion;
    }

    public double[] getOverdispersions() {
        return new double[] {this.quantifiedTretIPOverdispersion, this.quantifiedCtrlIPOverdispersion,
                             this.quantifiedTretINPUTOverdispersion, this.quantifiedCtrlINPUTOverdispersion};
    }

    public void setQuantifiedExpansionEffect(double quantifiedExpansionEffect) {
        this.quantifiedExpansionEffect = quantifiedExpansionEffect;
    }

    public double getQuantifiedExpansionEffect() {
        return this.quantifiedExpansionEffect;
    }
}
