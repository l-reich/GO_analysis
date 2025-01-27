import java.util.Objects;

public class IdPair {
    private final String first;
    private final String second;

    public IdPair(String first, String second) {
        // Ensure the pair is always stored in a consistent order
        if (first.compareTo(second) < 0) {
            this.first = first;
            this.second = second;
        } else {
            this.first = second;
            this.second = first;
        }
    }

    public String getFirst() {
        return first;
    }

    public String getSecond() {
        return second;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        IdPair idPair = (IdPair) o;
        return first.equals(idPair.first) && second.equals(idPair.second);
    }

    @Override
    public int hashCode() {
        return Objects.hash(first, second);
    }

    @Override
    public String toString() {
        return "(" + first + "," + second + ")";
    }
}