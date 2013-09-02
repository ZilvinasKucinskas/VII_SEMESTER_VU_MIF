import java.util.Random;

/**
 * @author Žilvinas Kučinskas
 *
 * Aprasymas:
 *
 * Si klase yra duomenu kanalas, i ji ieina tam tikri duomenys,
 *   yra klaidos tikimybe p, su kuria duomenys gali issikraipyti
 */

public class Channel {
    
    private int probability;
    private Random random;

    /*
     * Konstruktorius, iskart nustato iskraipymo(klaida) tikimybe
     *   0-100
     */
    public Channel(int probability) {
        this.probability = probability;
        this.random = new Random();
    }

    /*
     * Metodas siuncia pranesima kanalu su tikimybe iskraipyti signala
     *   taigi paimam matrica ir kiekvienam tos matricos elementui
     *   ziurime ar reikia iskraipyti ar ne, ir iskraipome.
     * Taip mes simuliuojame tikra gyvenimiska kanala, pvz Wireless
     */
    public int[][] send(int[][] message) {
        int messageRows = message.length,
            messageColumns = message[0].length;
        
        int[][] result = new int[messageRows][messageColumns];

        for (int i = 0; i < messageRows; i++) { // iteruojam per eilutes
            for (int j =0; j < messageColumns; j++) { // iteruojam per stulpelius
                if (this.isDistorted()) { //jei reikia iskraipyti
                    // pridedam vieneta ir mod 2 nes abecele binarine
                    result[i][j] = (message[i][j] + 1) % 2;
                } else {
                    // jei nereik iskraipyt, tai praeina toks pat bitukas
                    result[i][j] = message[i][j];
                }
            }
        }

        return result;
    }

    /*
     * Metodas grazina true jei ivyksta iskraipymas
     *   ir false jei neivyksta iskraipymas
     */
    private boolean isDistorted() {
        if (random.nextInt(100) < this.probability) {
            return true;
        } else {
            return false;
        }
    }
}