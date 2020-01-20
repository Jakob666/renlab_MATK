package DifferentialMethylation;

@FunctionalInterface
public interface FileRecordTask {
    Runnable createTask(double tretMeth, double ctrlMeth, double tretExp, double ctrlExp, int peakIdx);
}
