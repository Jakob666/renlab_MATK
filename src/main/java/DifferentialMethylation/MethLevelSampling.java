package DifferentialMethylation;

@FunctionalInterface
public interface MethLevelSampling {
    Runnable createTask(double prevValue, double curValue, int peakIdx);
}
