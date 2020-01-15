package DifferentialMethylation;

public class Parameters {
    private double tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                   tretMethylation, ctrlMethylation, tretBkgExp, tretNonPeakBkgExp, ctrlBkgExp, ctrlNonPeakBkgExp, nonspecificEnrichment;

    private double[] tretIPSizeFactors, tretINPUTSizeFactors, tretIPNonPeakSizeFactors, tretINPUTNonPeakSizeFactors,
                     ctrlIPSizeFactors, ctrlINPUTSizeFactors, ctrlIPNonPeakSizeFactors, ctrlINPUTNonPeakSizeFactors;

    public Parameters() {}

    public double getTretIPOverdispersion() {
        return this.tretIPOverdispersion;
    }

    public void setTretIPOverdispersion(double value) {
        this.tretIPOverdispersion = value;
    }

    public double getTretINPUTOverdispersion() {
        return this.tretINPUTOverdispersion;
    }

    public void setTretINPUTOverdispersion(double value) {
        this.tretINPUTOverdispersion = value;
    }

    public double getCtrlIPOverdispersion() {
        return this.ctrlIPOverdispersion;
    }

    public void setCtrlIPOverdispersion(double value) {
        this.ctrlIPOverdispersion = value;
    }

    public double getCtrlINPUTOverdispersion() {
        return this.ctrlINPUTOverdispersion;
    }

    public void setCtrlINPUTOverdispersion(double value) {
        this.ctrlINPUTOverdispersion = value;
    }

    public double getTretMethylation() {
        return this.tretMethylation;
    }

    public void setTretMethylation(double value) {
        this.tretMethylation = value;
    }

    public double getCtrlMethylation() {
        return this.ctrlMethylation;
    }

    public void setCtrlMethylation(double value) {
        this.ctrlMethylation = value;
    }

    public double getTretBkgExp() {
        return this.tretBkgExp;
    }

    public void setTretBkgExp(double value) {
        this.tretBkgExp = value;
    }

    public double getTretNonPeakBkgExp() {
        return this.tretNonPeakBkgExp;
    }

    public void setTretNonPeakBkgExp(double value) {
        this.tretNonPeakBkgExp = value;
    }

    public double getCtrlBkgExp() {
        return this.ctrlBkgExp;
    }

    public void setCtrlBkgExp(double value) {
        this.ctrlBkgExp = value;
    }

    public double getCtrlNonPeakBkgExp() {
        return this.ctrlNonPeakBkgExp;
    }

    public void setCtrlNonPeakBkgExp(double value) {
        this.ctrlNonPeakBkgExp = value;
    }

    public void setNonspecificEnrichment(double value) {
        this.nonspecificEnrichment = value;
    }

    public double getNonspecificEnrichment() {
        return this.nonspecificEnrichment;
    }

    public void setTretINPUTSizeFactors(double[] tretINPUTSizeFactors) {
        this.tretINPUTSizeFactors = tretINPUTSizeFactors;
    }

    public double[] getTretINPUTSizeFactors() {
        return this.tretINPUTSizeFactors;
    }

    public void setTretINPUTNonPeakSizeFactors(double[] tretINPUTNonPeakSizeFactors) {
        this.tretINPUTNonPeakSizeFactors = tretINPUTNonPeakSizeFactors;
    }

    public double[] getTretINPUTNonPeakSizeFactors() {
        return this.tretINPUTNonPeakSizeFactors;
    }

    public void setTretIPSizeFactors(double[] tretIPSizeFactors) {
        this.tretIPSizeFactors = tretIPSizeFactors;
    }

    public double[] getTretIPSizeFactors() {
        return this.tretIPSizeFactors;
    }

    public void setTretIPNonPeakSizeFactors(double[] tretIPNonPeakSizeFactors) {
        this.tretIPNonPeakSizeFactors = tretIPNonPeakSizeFactors;
    }

    public double[] getTretIPNonPeakSizeFactors() {
        return this.tretIPNonPeakSizeFactors;
    }

    public void setCtrlINPUTSizeFactors(double[] ctrlINPUTSizeFactors) {
        this.ctrlINPUTSizeFactors = ctrlINPUTSizeFactors;
    }

    public double[] getCtrlINPUTSizeFactors() {
        return this.ctrlINPUTSizeFactors;
    }

    public void setCtrlINPUTNonPeakSizeFactors(double[] ctrlINPUTNonPeakSizeFactors) {
        this.ctrlINPUTNonPeakSizeFactors = ctrlINPUTNonPeakSizeFactors;
    }

    public double[] getCtrlINPUTNonPeakSizeFactors() {
        return this.ctrlINPUTNonPeakSizeFactors;
    }

    public void setCtrlIPSizeFactors(double[] ctrlIPSizeFactors) {
        this.ctrlIPSizeFactors = ctrlIPSizeFactors;
    }

    public double[] getCtrlIPSizeFactors() {
        return this.ctrlIPSizeFactors;
    }

    public void setCtrlIPNonPeakSizeFactors(double[] ctrlIPNonPeakSizeFactors) {
        this.ctrlIPNonPeakSizeFactors = ctrlIPNonPeakSizeFactors;
    }

    public double[] getCtrlIPNonPeakSizeFactors() {
        return this.ctrlIPNonPeakSizeFactors;
    }
}
