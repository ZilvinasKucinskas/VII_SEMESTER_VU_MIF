import java.io.File;
import java.util.Scanner;

/**
 * @author Žilvinas Kučinskas
 *
 * Pagrindine klase, atlieka meniu funkcija ir apima kodavima,
 * siuntima kanalu, dekodavima aukstu lygiu. (Zemu lygiu suprogramuota atskirose
 * klasese)
 */
public class Main {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        String txtMessage = "message_before_encoding.txt";
        File fileMessage = new File(txtMessage);
        String txtEncoded = "message_encoded.txt";
        File fileEncoded = new File(txtEncoded);
        String txtSentThroughChannel = "message_sent.txt";
        File fileSent = new File(txtSentThroughChannel);
        String txtDecoded = "message_decoded.txt";
        File fileDecoded = new File(txtDecoded);
        int probability = 5;
        BinaryMatrixReader vectorReader;
        BinaryMatrixWriter vectorWriter;
        Channel channel;

        int[][] message, messageEncoded,messageReceived, messageDecoded;

        String line;
        int choice = -1;

        while (choice != 5) {
            System.out.println("Meniu:");
            System.out.println("0. Kanalo pralaidumo nustatymas");
            System.out.println("1. Kodavimas");
            System.out.println("2. Siuntimas kanalu");
            System.out.println("3. Dekodavimas");
            System.out.println("4. HELP");
            System.out.println("5. EXIT");
            System.out.println();

            try {
                line = scanner.nextLine();
                choice = Integer.parseInt(line);

                switch (choice) {
                    case 0:
                        System.out.println("PASIRINKTAS KODAVIMAS");
                        System.out.println("Prasome ivesti klaidos tikimybe p (0-100)");
                        line = scanner.nextLine();
                        probability = Integer.parseInt(line);
						System.out.println("Kanalo pralaidumas nustatytas");
                        break;
                    case 1:
                        System.out.println("PASIRINKTAS KODAVIMAS");
                        System.out.println("Skaitoma is failo " + txtMessage);
                        vectorReader = new BinaryMatrixReader(fileMessage);
                        message = vectorReader.readVector();
                        if (message != null) {
                            System.out.println("Nuskaiteme faila, koduojame");
                            messageEncoded = Encoder.encode(message);
                            vectorWriter = new BinaryMatrixWriter(fileEncoded);
                            boolean success = vectorWriter.writeVector(messageEncoded);
                            if (success) {
                                System.out.println("Pavyko uzkoduoti zinute");
                                System.out.println("Irasyta i faila " + txtEncoded);
                            } else {
                                System.out.println("Nepavyko irasyti i faila");
                            }
                        } else {
                            System.out.println("Nenuskaiteme failo");
                        }
                        break;
                    case 2:
                        System.out.println("PASIRINKTAS SIUNTIMAS KANALU");
                        System.out.println("Skaitoma is failo " + txtEncoded);
                        vectorReader = new BinaryMatrixReader(fileEncoded);
                        messageEncoded = vectorReader.readVector();
                        if (messageEncoded != null) {
                            System.out.println("Nuskaiteme faila, siunciame kanalu");
                            channel = new Channel(probability);
                            messageReceived = channel.send(messageEncoded);
                            vectorWriter = new BinaryMatrixWriter(fileSent);
                            boolean success = vectorWriter.writeVector(messageReceived);
                            if (success) {
                                System.out.println("Pavyko persiusti zinute kanalu");
                                System.out.println("Irasyta i faila " + txtSentThroughChannel);
                                System.out.println("Faila galima atsidarius txt faila modifikuoti pries dekodavima");
                            } else {
                                System.out.println("Nepavyko irasyti i faila");
                            }
                        } else {
                            System.out.println("Nenuskaiteme failo");
                        }
                        break;
                    case 3:
                        System.out.println("PASIRINKTAS DEKODAVIMAS");
                        System.out.println("Skaitoma is failo " + txtSentThroughChannel);
                        vectorReader = new BinaryMatrixReader(fileSent);
                        messageReceived = vectorReader.readVector();
                        if (messageReceived != null) {
                            System.out.println("Nuskaiteme faila, dekoduojame");
                            messageDecoded = Decoder.decode(messageReceived);
                            if (messageDecoded == null) {
                                System.out.println("Nepavyko dekoduoti, prasau "
                                        + "persiusti is naujo kanalu zinute");
                                break;
                            }
                            vectorWriter = new BinaryMatrixWriter(fileDecoded);
                            boolean success = vectorWriter.writeVector(messageDecoded);
                            if (success) {
                                System.out.println("Pavyko dekoduoti zinute");
                                System.out.println("Irasyta i faila " + txtDecoded);
                            } else {
                                System.out.println("Nepavyko irasyti i faila");
                            }
                        } else {
                            System.out.println("Nenuskaiteme failo");
                        }
                        break;
                    case 4:
                        System.out.println("PASIRINKTAS HELP");
                        System.out.println("Standartiniai parametrai:");
                        System.out.println("Klaidos tikimybe "
                                + probability);
                        System.out.println("Pradinis vektorius saugomas faile "
                                + txtMessage);
                        System.out.println("Uzkoduotas vektorius saugomas faile "
                                + txtEncoded);
                        System.out.println("Persiustas kanalu vektorius saugomas faile "
                                + txtSentThroughChannel);
                        System.out.println("Dekoduotas vektorius saugomas faile "
                                + txtDecoded);
                        break;
                    case 5:
                        System.out.println("Programa baigta");
                        break;
                    default:
                        System.out.println("Tikriausiai ivedimo klaida");

                        }
            } catch (Exception e) {
                System.out.println("Ivyko klaida" + e.getMessage());
                e.printStackTrace();
            }
        }
    }
}