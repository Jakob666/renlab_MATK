package DifferentialMethylation;

public class Parameters {
    private double tretIPOverdispersion, tretINPUTOverdispersion, ctrlIPOverdispersion, ctrlINPUTOverdispersion,
                   tretMethylation, ctrlMethylation, tretBkgExp, tretNonPeakBkgExp, ctrlBkgExp, ctrlNonPeakBkgExp, nonspecificEnrichment;

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
}
