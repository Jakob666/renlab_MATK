package DifferentialMethylation;

@FunctionalInterface
public interface CreateTask {
    Runnable getTask(int geneIdx);
}
