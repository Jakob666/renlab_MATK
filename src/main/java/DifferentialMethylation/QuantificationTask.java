package DifferentialMethylation;

@FunctionalInterface
public interface QuantificationTask {
    Runnable createTask(int peakIdx);
}
