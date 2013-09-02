/**
 * @author Žilvinas Kučinskas
 * 
 * Aprasymas:
 *
 * Klase skirta aprasyti operacijas atliekamas su matricomis,
 *   tokias kaip daugyba, sudetis
 */

public class BinaryMatrixOperations {

    /*
     * Metodas sudaugina 2 matricas
     *   [m x n] x [n x k] = [m, k]
     * Remiamasi, jog abecele yra {0, 1}
     */
    public static int[][] multiply(int a[][], int b[][]) {
        int aRows = a.length,
            aColumns = a[0].length,
            bRows = b.length,
            bColumns = b[0].length;

        // A matricos stulpeliu turi but tiek kiek B 
        if (aColumns != bRows) {
            throw new IllegalArgumentException("A matricos eiluciu: "
                    + aColumns + " neatitinka B matricos stulpeliu "
                    + bRows + ".");
        }

        int[][] result = new int[aRows][bColumns];

        for (int i = 0; i < aRows; i++) { // iteruojam per A eilutes
            for (int j = 0; j < bColumns; j++) { // iteruojam per B stulpelius
                for (int k = 0; k < aColumns; k++) { // iteruojam per A stulpelius
                    result[i][j] += a[i][k] * b[k][j];
                }
                result[i][j] %= 2; // operuojam abecele {0, 1}
            }
        }

        return result;
    }


    /*
     * Metodas sudeda dvi matricas
     *    [m x n] + [m x n] -> dimensijos turi buti vienodos
     * Remiamasi jog abecele yra binarine [0, 1]
     */
    public static int[][] addition(int[][] a, int[][] b) {
        int aRows = a.length,
            aColumns = a[0].length,
            bRows = b.length,
            bColumns = b[0].length;

        if (aRows != bRows || aColumns != bColumns) {
            throw new IllegalArgumentException("Matricos A eiluciu skaicius "
                    + aRows +  " neatitinka B eiluciu skaiciaus " + bRows
                    + " arba matricos A stulpeliu skaicius " + aColumns
                    + " neatitinka B stulpeliu skaiciaus " + bColumns + ".");
        }

        int[][] result = new int[aRows][aColumns];

        for (int i = 0; i < aRows; i++) { // iteruojam per eilutes
            for (int j = 0; j < aColumns; j++) { // iteruojam per stulpelius
                // skaiciuojam rezultata ir dar darom mod 2
                // nes abecele binarine
                result[i][j] = (a[i][j] + b[i][j]) % 2;
            }
        }

        return result;
    }
}