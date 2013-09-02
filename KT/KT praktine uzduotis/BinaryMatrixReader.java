import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

/**
 * @author Žilvinas Kučinskas
 *
 * Aprasymas:
 *
 * Klase skirta skaityti binarini vektoriu is failo kaip stringa
 *
 */

public class BinaryMatrixReader {
    private File file;

    /*
     * Konstruktorius, kuris kaip parametra gauna faila is kurio skaitysime
     */
    public BinaryMatrixReader(File file) {
        this.file = file;
    }

    /*
     * Metodas skirtas nuskaityti vektoriu is failo,
     *   kuris tame faile yra Stringas, automatiskai pavercia i matrica [1xn],
     *   kur n yra vektoriaus ilgis
     */
    public int[][] readVector(){
        int vector[][] = null;
        
        try {
            FileReader fstream = new FileReader(file);
            BufferedReader in = new BufferedReader(fstream);
            String vectorString = in.readLine();

            vectorString.replaceAll("//s", ""); //remove all whitespaces
            vector = new int[1][vectorString.length()];

            for (int i = 0; i < vectorString.length(); i++) {
                int value = Character.getNumericValue(vectorString.charAt(i));
                if (value != 0 && value != 1) {
                    throw new IllegalArgumentException("Rastas simbolis, kuris"
                            + " nera vienetas arba nulis");
                }
                vector[0][i] = value;
            }

            in.close();
        } catch (Exception e) {
            System.out.println("Error: " + e);
            return null;
        }
        return vector;
    }
}