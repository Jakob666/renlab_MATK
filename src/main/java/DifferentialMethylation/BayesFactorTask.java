package DifferentialMethylation;

@FunctionalInterface
public interface BayesFactorTask {
    Runnable createTask(int peakIdx);
}
