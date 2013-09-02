/**
 * @author Žilvinas Kučinskas
 *
 * Aprasymas:
 *
 * Golejaus kodo (23, 12, 7) dekodavimo veiksmai atliekami sioje klaseje
 */
public class Decoder {

    /*
     * Metodas naudojasi dekodavimo algoritmu 3.7.1
     *   Grazina dekoduota zinute jeigu yra 3 klaidos ir maziau,
     *   jei daugiau klaidu tai grazina null, ir tai reiskia,
     *   kad mum reikia is naujo kanalu siusti zinute
     */
    public static int[][] decode(int[][] message) {
        int[][] decodedMessageC23 = new int[1][12]; // cia jau grazinsim dekoduota zinute
        // klaidu patternas, jis bus 1 eil ir 24 stulp golejaus kode
        int[][] mistakesPattern = new int[1][24];
        // Prie gauto zodzio pridedam viena bituka, privedam prie formos
        // is C23 i C24
        int[][] messageC24 = new int[1][24];
        // kopijuojam turini zinutes
        System.arraycopy(message[0], 0, messageC24[0], 0, 23);
        // jei vienetuku lyginis skaicius pridedam nuliuka,
        // o jei nelyginis tai pridedam nuliuka kaip 24 simboli
        if (Decoder.weight(message) % 2 == 0) {
            messageC24[0][23] = 1;
        } else {
            messageC24[0][23] = 0;
        }

        /* TESTAVIMAS
        // pasiskaiciuojam sindroma, jis bus 12 ilgio
        int[][] u1 = new int[1][12];
        System.arraycopy(messageC24[0], 0, u1[0], 0, 12);
        int[][] u2 = new int[1][12];
        System.arraycopy(messageC24[0], 12, u2[0], 0, 12);
        //int[][] syndrome = BinaryMatrixOperations.addition(u1, Decoder.syndrome(u2));
        int[][] syndrome = Decoder.syndrome(u2);
         * *
         */

        //pasiskaiciuojam sindroma s=wH=u1+u2B ir u2*B iskelta i metoda syndrome
        int[][] u1 = new int[1][12];
        System.arraycopy(messageC24[0], 0, u1[0], 0, 12);
        int[][] u2 = new int[1][12];
        System.arraycopy(messageC24[0], 12, u2[0], 0, 12);
        int[][] syndrome = BinaryMatrixOperations.addition(u1, Decoder.syndrome(u2));

        /* TESTAVIMAS
        for (int i = 0; i < u1[0].length; i++) {
            System.out.print(u1[0][i]);
        }
        System.out.println();
        for (int i = 0; i < u2[0].length; i++) {
            System.out.print(u2[0][i]);
        }
        System.out.println();
        for (int i = 0; i < syndrome[0].length; i++) {
            System.out.print(syndrome[0][i]);
        }
        System.out.println();
        */
        
        // pirmas zingsnis
        // jei sindromo svoris <= 3
        if (Decoder.weight(syndrome) <= 3) {
            // tai u = [s, 0]
            System.arraycopy(syndrome[0], 0, mistakesPattern[0], 0, 12);
            for (int i = 12; i < 24; i++) {
                mistakesPattern[0][i] = 0;
            }
            for (int i = 0; i < mistakesPattern[0].length; i++) {
                if (mistakesPattern[0][i] == 1) {
                    System.out.println("Klaida ivyko pozicijoje " + i);
                }
            }
            // paskaiciuojam prideje klaidu patterna prie gauto zodzio
            int[][] decodedMessageC24 = BinaryMatrixOperations.addition(messageC24, mistakesPattern);
            // ir grazinam jau dekoduota rezultata, kuris yra 12 pirmuju bituku
            System.arraycopy(decodedMessageC24[0], 0, decodedMessageC23[0], 0, 12);
            return decodedMessageC23;
        }

        int[][] bi; // vienetines matricos eilute tam tikroj pozicijoj
        // antrasis zingsnis
        // kiekvienai eilutei is B
        for (int i = 0; i < 12; i++) {
            bi = new int[1][12];
            System.arraycopy(MatrixData.B[i], 0, bi[0], 0, 12);
            // jei svoris sindromo + bi eilutes <= 2
            if (Decoder.weight(BinaryMatrixOperations.addition(syndrome, bi)) <= 2) {
                // tai u = [s+bi, ei], cia ei tai vienetines matricos eilute
                //   pozicijoje i
                int[][] result = BinaryMatrixOperations.addition(syndrome, bi);
                System.arraycopy(result[0], 0, mistakesPattern[0], 0, 12);
                System.arraycopy(MatrixData.I[i], 0, mistakesPattern[0], 12, 12);
                for (int j = 0; j < mistakesPattern[0].length; j++) {
                    if (mistakesPattern[0][j] == 1) {
                        System.out.println("Klaida ivyko pozicijoje " + j);
                    }
                }
                // paskaiciuojam prideje klaidu patterna prie gauto zodzio
                int[][] decodedMessageC24 = BinaryMatrixOperations.addition(messageC24, mistakesPattern);
                // ir grazinam jau dekoduota rezultata, kuris yra 12 pirmuju bituku
                System.arraycopy(decodedMessageC24[0], 0, decodedMessageC23[0], 0, 12);
                return decodedMessageC23;
            }
        }

        // pasiskaiciuojam antraji sindroma sB, jis bus 12 ilgio
        int[][] syndrome2 = BinaryMatrixOperations.multiply(syndrome, MatrixData.B);

        // 3 zingsnis
        // jei paskaiciuoto sindromo svoris maziau uz 3 tai
        if (Decoder.weight(syndrome2) <= 3) {
            // u = [0, sB]
            for (int i = 0; i < 12; i++) {
                mistakesPattern[0][i] = 0;
            }
            System.arraycopy(syndrome2[0], 0, mistakesPattern[0], 12, 12);
            for (int i = 0; i < mistakesPattern[0].length; i++) {
                if (mistakesPattern[0][i] == 1) {
                    System.out.println("Klaida ivyko pozicijoje " + i);
                }
            }
            // paskaiciuojam prideje klaidu patterna prie gauto zodzio
            int[][] decodedMessageC24 = BinaryMatrixOperations.addition(messageC24, mistakesPattern);
            // ir grazinam jau dekoduota rezultata, kuris yra 12 pirmuju bituku
            System.arraycopy(decodedMessageC24[0], 0, decodedMessageC23[0], 0, 12);
            return decodedMessageC23;
        }

        // 4 zingsnis
        for (int i = 0; i < 12; i++) {
            bi = new int[1][12];
            System.arraycopy(MatrixData.B[i], 0, bi[0], 0, 12);
            // jei svoris sindromo sB + bi eilutes <= 2
            if (Decoder.weight(BinaryMatrixOperations.addition(syndrome2, bi)) <= 2) {
                // tai u = [ei, sB + bi], cia ei tai vienetines matricos eilute
                //   pozicijoje i
                int[][] result = BinaryMatrixOperations.addition(syndrome2, bi);
                System.arraycopy(MatrixData.I[i], 0, mistakesPattern[0], 0, 12);
                System.arraycopy(result[0], 0, mistakesPattern[0], 12, 12);
                for (int j = 0; j < mistakesPattern[0].length; j++) {
                    if (mistakesPattern[0][j] == 1) {
                        System.out.println("Klaida ivyko pozicijoje " + j);
                    }
                }
                // paskaiciuojam prideje klaidu patterna prie gauto zodzio
                int[][] decodedMessageC24 = BinaryMatrixOperations.addition(messageC24, mistakesPattern);
                // ir grazinam jau dekoduota rezultata, kuris yra 12 pirmuju bituku
                System.arraycopy(decodedMessageC24[0], 0, decodedMessageC23[0], 0, 12);
                return decodedMessageC23;
            }
        }
        // prasome pakartoti siuntima jei neradom klaidu patterno
        return null;
    }


    /*
     * Metodas skirtas apskaiciuoti sindroma
     *     sindromas = uzkoduoto žodžio 12bitu * B
     *      kiti reikalingi veiksmai atliekami decode metode
     */
    private static int[][] syndrome(int[][] message) {
        return BinaryMatrixOperations.multiply(message, MatrixData.B);
    }

    /*
     * Metodas skirtas paskaiciuoti zodzio svoriui
     */
    private static int weight(int[][] message) {
        int result = 0;

        for (int i = 0; i < message[0].length; i++) {
            if (message[0][i] != 0) {
                result++;
            }
        }

        return result;
    }
}
