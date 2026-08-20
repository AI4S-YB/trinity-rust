import jaligner.NeedlemanWunschGotohBanded;
import jaligner.NeedlemanWunschGotoh;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixGenerator;
public class JAlignGolden {
    public static void main(String[] args) throws Exception {
        Matrix m = MatrixGenerator.generate(4f, -5f);
        java.io.BufferedReader r = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
        String line;
        while ((line = r.readLine()) != null) {
            if (line.isEmpty()) continue;
            String[] p = line.split("\t");
            String mode = p[0]; String s1 = p[1]; String s2 = p[2];
            Sequence seq1 = new Sequence("A", s1);
            Sequence seq2 = new Sequence("B", s2);
            Alignment a;
            if (mode.equals("B")) { int bw = Integer.parseInt(p[3]); a = NeedlemanWunschGotohBanded.align(seq1, seq2, m, 10f, 1f, bw); }
            else { a = NeedlemanWunschGotoh.align(seq1, seq2, m, 10f, 1f); }
            System.out.println(String.valueOf(a.getSequence1()) + "\t" + String.valueOf(a.getSequence2()) + "\t" + a.getScore() + "\t" + a.getStart1() + "\t" + a.getStart2());
        }
    }
}
