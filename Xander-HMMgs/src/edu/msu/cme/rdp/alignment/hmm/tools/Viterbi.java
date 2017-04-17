/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.alignment.hmm.tools;

import edu.msu.cme.rdp.alignment.hmm.scoring.ViterbiScorer;
import edu.msu.cme.rdp.alignment.hmm.scoring.HMMScorer;
import edu.msu.cme.rdp.alignment.hmm.HMMER3bParser;
import java.io.File;
import java.util.Arrays;
import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.alignment.hmm.XSC;
import edu.msu.cme.rdp.alignment.hmm.XSTATES;
import edu.msu.cme.rdp.readseq.readers.Sequence;

import static edu.msu.cme.rdp.alignment.hmm.TSC.*;
import static edu.msu.cme.rdp.alignment.hmm.XSTATES.*;

import static java.lang.StrictMath.*;

/**
 *
 * @author fishjord
 */
public class Viterbi {

    private static final int M = 0, I = 1, D = 2;

    public static double viterbi(Sequence seq, ProfileHMM hmm) {
        char[] bases = seq.getSeqString().toCharArray();

        double[][][] matrix = new double[3][bases.length + 1][hmm.M() + 1];

        double[][] xmx = new double[bases.length + 1][XSTATES.values().length];

        xmx[0][N.ordinal()] = 0;
        xmx[0][B.ordinal()] = hmm.xsc(N, XSC.MOVE);
        xmx[0][XSTATES.E.ordinal()] = xmx[0][C.ordinal()] = xmx[0][J.ordinal()] = Double.NEGATIVE_INFINITY;

        for (int k = 0; k <= hmm.M(); k++) {
            matrix[M][0][k] = matrix[I][0][k] = matrix[D][0][k] = Double.NEGATIVE_INFINITY;
        }

        double esc = Double.NEGATIVE_INFINITY;

        for (int i = 1; i < bases.length + 1; i++) {

//            if(i > 5) System.exit(1);
            char xi = bases[i - 1];
            double sc;

            matrix[M][i][0] = matrix[I][i][0] = matrix[D][i][0] = Double.NEGATIVE_INFINITY;
            xmx[i][XSTATES.E.ordinal()] = Double.NEGATIVE_INFINITY;

            for (int k = 1; k <= hmm.M(); k++) {
                sc = max(matrix[M][i - 1][k - 1] + hmm.tsc(k - 1, MM),
                        matrix[I][i - 1][k - 1] + hmm.tsc(k - 1, IM));
                sc = max(sc, matrix[D][i - 1][k - 1] + hmm.tsc(k - 1, DM));
                //System.out.println("SC = " + sc + " XMX: " + xmx[i - 1][B.ordinal()] + " TSC[" + (k - 1) + "][BM]= " + hmm.tsc(k - 1, BM));
                sc = max(sc, xmx[i - 1][B.ordinal()] + hmm.tsc(k - 1, BM));
                //System.out.println("SC = " + sc + " XMX: " + xmx[i - 1][B.ordinal()] + " TSC[" + (k - 1) + "][BM]= " + hmm.tsc(k - 1, BM));

                matrix[M][i][k] = sc + hmm.msc(k, xi);

                //xmx[i][XSTATES.E.ordinal()] = max(xmx[i][XSTATES.E.ordinal()], matrix[M][i][k] + esc);

                sc = max(matrix[M][i - 1][k] + hmm.tsc(k, MI),
                        matrix[I][i - 1][k] + hmm.tsc(k, II));
                matrix[I][i][k] = sc + hmm.isc(k, xi);

                matrix[D][i][k] = max(matrix[M][i][k - 1] + hmm.tsc(k - 1, MD),
                        matrix[D][i][k - 1] + hmm.tsc(k - 1, DD));
            }
            sc = max(matrix[M][i - 1][hmm.M() - 1] + hmm.tsc(hmm.M() - 1, MM),
                    matrix[I][i - 1][hmm.M() - 1] + hmm.tsc(hmm.M() - 1, IM));
            sc = max(sc, matrix[D][i - 1][hmm.M() - 1] + hmm.tsc(hmm.M() - 1, DM));
            sc = max(sc, xmx[i - 1][B.ordinal()] + hmm.tsc(hmm.M() - 1, BM));

            matrix[M][i][hmm.M()] = sc + hmm.msc(hmm.M(), xi);

            matrix[D][i][hmm.M()] = max(matrix[M][i][hmm.M() - 1] + hmm.tsc(hmm.M() - 1, MD),
                    matrix[D][i][hmm.M() - 1] + hmm.tsc(hmm.M() - 1, DD));

            //E
            sc = max(xmx[i][XSTATES.E.ordinal()], matrix[M][i][hmm.M()]);
            xmx[i][XSTATES.E.ordinal()] = max(sc, matrix[D][i][hmm.M()]);

            //J
            sc = xmx[i - 1][J.ordinal()] + hmm.xsc(J, XSC.LOOP);
            xmx[i][J.ordinal()] = max(sc, xmx[i][XSTATES.E.ordinal()] + hmm.xsc(XSTATES.E, XSC.LOOP));

            //C
            sc = xmx[i - 1][C.ordinal()] + hmm.xsc(C, XSC.LOOP);
            xmx[i][C.ordinal()] = max(sc, xmx[i][XSTATES.E.ordinal()] + hmm.xsc(XSTATES.E, XSC.MOVE));

            //N
            xmx[i][N.ordinal()] = xmx[i - 1][N.ordinal()] + hmm.xsc(N, XSC.LOOP);

            //B
            sc = xmx[i][N.ordinal()] + hmm.xsc(N, XSC.MOVE);
            xmx[i][B.ordinal()] = max(sc, xmx[i][J.ordinal()] + hmm.xsc(J, XSC.MOVE));

            /*System.out.print("xmx[" + i + "][E]= " + xmx[i][XSTATES.E.ordinal()] + " ");
            System.out.print("xmx[" + i + "][J]= " + xmx[i][XSTATES.J.ordinal()] + " ");
            System.out.print("xmx[" + i + "][C]= " + xmx[i][XSTATES.C.ordinal()] + " ");
            System.out.print("xmx[" + i + "][N]= " + xmx[i][XSTATES.N.ordinal()] + " ");
            System.out.println("xmx[" + i + "][B]= " + xmx[i][XSTATES.B.ordinal()] + " ");*/
        }


        System.out.println("MMX[" + bases.length + "][" + hmm.M() + "]= " + matrix[M][bases.length][hmm.M()]);
        System.out.println("IMX[" + bases.length + "][" + hmm.M() + "]= " + matrix[I][bases.length][hmm.M()]);
        System.out.println("DMX[" + bases.length + "][" + hmm.M() + "]= " + matrix[D][bases.length][hmm.M()]);
        System.out.println("XMX[L][C]: " + (xmx[bases.length][C.ordinal()] + hmm.xsc(C, XSC.MOVE)));
        System.out.println("XMX[L][E]: " + (xmx[bases.length][XSTATES.E.ordinal()]));

        return xmx[bases.length][C.ordinal()] + hmm.xsc(C, XSC.MOVE);
    }

    public static void main(String[] args) throws Exception {
        String testModel = "/work/fishjord/other_projects/hmm_graph/testing_models/nifH/nifH.hmm";
        //String testModel = "/work/fishjord/other_projects/hmm_graph/testing_models/nifH/reversed_nifH.hmm";
        //String testModel = "/work/fishjord/other_projects/hmm_graph/testing_models/tiny_nifh/tiny_nifh.hmm";

        //Sequence testSeq = new Sequence("AAB63257", "", "mrqvaiygkggigkstttqnltaglgemgkkimivgcdpkadstrlvlgglaqktvldtlreegedieldtvlkvgyagikgvesggpepasaagrgiitsigllerlgayeadldyvfydvlgdvvcggfampiregkaqeiyivcsaemmglyaanniakgiskyantggvrlgglicnsrkvdgeadlvsrvakeigtqmihfvpatmrcrrrksikrqlstfrpmtqadeyrtlarkidgndmfvvprpmsidrleailmehgild");
        //Sequence testSeq = new Sequence("pass20_aligned-sample1", "", "PLWLASANTMPVGTQEGWQAYLSQSDKAMRQVVEGRPISGKGVIGQETAEEIVLIVTITNTGKRVGPKGAQMKADSLRVIANGLFQDTILDVAGILGKSGNVEVDDSSTEKQTGISSAELDSPCPGIGCPFHGVLTADESVEEKRVGADIDGVMLHLLADPVCTGFEQPIRIFLAQRVVIIMTGEMWAVYKAGNIAENIVSFAHLTGIRLGGIITNKMPVRARELGICDRTSENAGTKMNFDIIRDVYIQAKEKRSMDVPSKEWKEEVAAEDARSVRHLTKMETGVFSFPKPHDVDSYVLVDAGEMRAEDNSLFSHADLNSSHGS");
        //Sequence testSeq = new Sequence("AAB63257", "", "mrqvaiygkggigkstttqnltaglgemgkkimivgcdpkadstrlvlgglaqktvldtlreegedieldtvlkvgyagikgvesggpepasaagrgiitsigllerlgayeadldyvfydvlgdvvcggfampiregkaqeiyivcsaemmglyaanniakgiskyantggvrlgglicnsrkvdgeadlvsrvakeigtqmihfvpatmrcrrrksikrqlstfrpmtqa");
        //Sequence testSeq = new Sequence("test", "", "MAGLRQIAFYGKGGIGKSTTSQNTLAALVDLGQKILIVGCDPKADSTRLILNAKAQDTVLHLAAKEGSVEDLEVEDVLKVGYKGIKCVESGGPEPGVGCAGRGVITSINFLEENGAYDDVDYVSYDVLGDVVCGGFAMPIRENKAQEIYIVMSGEMMALYAANNIAKGILKYAHSGGVRLGGLICNERQTDRELDLAEALAAKLNSRLIHFVPRDNIVQHAELRKMTVIQYAPESQQAAEYRALADKIHANSGQGTVPTPITMEELEDMLLDFGVMKTDEQMLAELQAKEAAAAAQ");
        //Sequence testSeq = new Sequence("test", "", "nvwdfyhagithasaqwtgvdeppgagllrnrynqdyltllgdyghaisgplaktaftasrpafdnawrdtpealkasgrsvtaagnvhvfpnmwiqqnfwmvalrmpmgpmkteiwwytfvteemdeetsdhivkrssrhngpagmtelddgen");

        //Sequence testSeq = new Sequence("", "", "EAKLRQIAIYGKGGIGKSTTSQNTLAALAELGKKVLIVGCDPKADSTRLILHAKAQDTVLDLAAEKGSVEDLELEDVLKEGYKDIKCVESGGPEPGVGCAGRGVITSINFLEEEGAYEDLDYVSYDVLGDVVCGGFAMPIRENKAQEIYIVVSGEMMALYAANNIAKGILKYAKSGGVRLGGLICNSRKVDREKELIEELAKKLGTKLIHFVPRDNIVQKAELRRKTVIEYAPESKQADEYRELAKKIVENEKLVIPTPITMDELEELLLEFGILEEEDESLVEKAAAEAAA");
        Sequence testSeq = new Sequence("", "", "MMKMRQIAFYGKGGIGKSTTTQNTMAAMAEMGQKIMIVGCDPKADSTRLILHSKAQDTVMHMAAEKGSVEDLELDDVMKKGYKDIKCVESGGPEPGVGCAGRGVITSINFLEENGAYEDVDYVFYDVLGDVVCGGFAMPIRENKAQEIYIVMSGEMMAMYAANNICKGILKYAHSGGVRLGGLICNSRKTDREQELIEELAKKLNTQMIHFVPRDNIVQHAEIRRMTVIEYDPDCKQADEYRQLAKKIHNNDMKVIPTPITMDELEDMLMEFGIMKQEDESIVEKTAAEAKA");
        //Sequence testSeq = new Sequence("", "", "pagraarparrvkkqkkviywtqllksetvskqllttccgitgvlctnnkvayngrekryktirvnlisntykskdyrkfeticenkkigiidcekvyifkcfyhyrkinvykflneymevyylgattqynasnffeqktwtcnkknkvekirvhgmqclscpfdinnyyvikeeyisiiklilnknryrrntiktyskyvinysarlecqismlknrrrdnyccykkcnksvindynkkifyglnlkhtyinmktlkvirkkvlilkpniywinlskintfltyvnesrrvekvrlyqysdiseslryykqlkklkitnrlyktkkrfkmniyvkrinneyrilkwylsnikstliffvkrqmvnlriifdevlkqrknearyygnyykeeenarmqsiciirerkiwnrpykiyeqqiliksgkgynygyyiiytvhlyrphrsryiqkntlkmtyisilccfnciptvpsle");
        //Sequence testSeq = new Sequence("", "", "eengayedidyvsydvlgdvvcggfampirenkaqeiyivcsgemmamyaanniakgivkyansggvrlaglicnsrntdredeliealasklgtqmihfvprdnavqhaeirrmtvieydpkhkqadeyrqlamkivnntkfviptpiemeeleellmefgimevedes");
        //Sequence testSeq = new Sequence("", "", new StringBuilder("eengayedidyvsydvlgdvvcggfampirenkaqeiyivcsgemmamyaanniakgivkyansggvrlaglicnsrntdredeliealasklgtqmihfvprdnavqhaeirrmtvieydpkhkqadeyrqlamkivnntkfviptpiemeeleellmefgimevedes").reverse().toString());
        
        ProfileHMM model = HMMER3bParser.readUnnormalized(new File(testModel));
        //ProfileHMM model = HMMER3bParser.readModel(new File(testModel));
        double expectedScore = 421.4f;

        System.out.println(testSeq.getSeqString().length());
        System.out.println(model.M());

        model.reconfigureLength(testSeq.getSeqString().length());

        int L = testSeq.getSeqString().length();
        double p1 = (double) L / (L + 1);
        double null1 = L * log(p1) + log(1 - p1);
        System.out.println("null1: " + null1);

        double d = Viterbi.viterbi(testSeq, model);
        System.out.println(d);
        System.out.println((d - null1) / log(2));

        System.out.println();

        HMMScorer scorer = new ViterbiScorer(model, 1);

        double maxScore = Double.NEGATIVE_INFINITY;
        int bestPos = 0;
        for (char c : testSeq.getSeqString().toCharArray()) {
            scorer.consume(c);
            double s = scorer.getScore();
            
            System.out.println(s + /*"\t" + scorer.getBestPossibleScore() +*/ "\t" + Arrays.asList(scorer.getMaxScore()));

            if (s > maxScore) {
                maxScore = s;
                bestPos = scorer.getQueryLength();
            }
        }

        System.out.println(scorer.getScore());

        System.out.println("Best score= " + maxScore + " bestPos= " + bestPos);
        System.out.println("Score= " + scorer.getScore());
    }
}
