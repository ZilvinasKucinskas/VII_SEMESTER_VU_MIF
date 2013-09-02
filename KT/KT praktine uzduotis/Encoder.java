/**
 * @author Žilvinas Kučinskas
 *
 * Aprasymas:
 *
 * Golejaus kodo (23, 12, 7) kodavimo veiksmai atliekami sioje klaseje
 */

public class Encoder {
    
    /*
     * Metodas uzkoduoti vektoriu golejaus kodas
     */
    public static int[][] encode(int[][] message) {
        return BinaryMatrixOperations.multiply(message, MatrixData.G);
    }
}