package edu.cmu.sphinx.tools.corpusBrowser;

import edu.cmu.sphinx.tools.batch.BatchNISTRecognizer;
import edu.cmu.sphinx.tools.corpus.Corpus;
import edu.cmu.sphinx.tools.corpus.Utterance;
import edu.cmu.sphinx.tools.corpus.Word;
import edu.cmu.sphinx.tools.corpus.xml.CorpusXMLReader;
import edu.cmu.sphinx.util.props.ConfigurationManager;
import edu.cmu.sphinx.util.props.PropertyException;

import javax.swing.*;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeSelectionModel;
import java.awt.*;
import java.io.*;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.TreeSet;

import com.intellij.uiDesigner.core.GridLayoutManager;
import com.intellij.uiDesigner.core.GridConstraints;
import com.intellij.uiDesigner.core.Spacer;

/**
 * Copyright 1999-2006 Carnegie Mellon University.
 * Portions Copyright 2002 Sun Microsystems, Inc.
 * Portions Copyright 2002 Mitsubishi Electric Research Laboratories.
 * All Rights Reserved.  Use is subject to license terms.
 * <p/>
 * See the file "license.terms" for information on usage and
 * redistribution of this file, and for a DISCLAIMER OF ALL
 * WARRANTIES.
 * <p/>
 * User: Peter Wolf
 * Date: Feb 2, 2006
 * Time: 9:40:58 AM
 */
public class CorpusBrowser implements TreeSelectionListener {

    private ConfigurationManager cm;
    private Corpus corpus;

    private JTree phonemeTree;
    private JList utterances;
    private JList words;
    private JList characters;
    private JTree characterTree;
    private JTree wordTree;
    private JTree utteranceTree;
    private JSlider offset;
    private JSlider bandwidth;
    private JScrollPane phonemeScroll;
    private JScrollPane characterScroll;
    private JScrollPane wordScroll;
    private JScrollPane utteranceScroll;
    private JPanel mainPanel;
    private Font font;

    CorpusBrowser(String fontFile) {
        try {
            font = Font.createFont(Font.TRUETYPE_FONT, new FileInputStream(fontFile)).deriveFont(20.0f);
        } catch (FontFormatException e) {
            throw new Error(e);
        } catch (FileNotFoundException e) {
            throw new Error(e);
        } catch (IOException e) {
            throw new Error(e);
        }
    }

    public static void main(String[] args) {

        if (args.length != 3) {
            System.out.println(
                    "Usage: CorpusBrowser propertiesFile corpusFile fontFile");
            System.exit(1);
        }

        String propertiesFile = args[0];
        String corpusFile = args[1];
        String fontFile = args[2];

        try {
            URL url = new File(propertiesFile).toURI().toURL();
            ConfigurationManager cm = new ConfigurationManager(url);

            Corpus corpus = new CorpusXMLReader( new FileInputStream(corpusFile) ).read();

            final CorpusBrowser cb = new CorpusBrowser(fontFile);
            JFrame f = new JFrame("CorpusBrowser");
            cb.setCorpus(cm, corpus);

            f.setContentPane(cb.mainPanel);
            f.pack();
            f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

            f.setVisible(true);


        } catch (MalformedURLException e) {
            throw new Error(e);
        } catch (PropertyException e) {
            throw new Error(e);
        } catch (IOException e) {
            throw new Error(e);
        }
    }


    private void setCorpus(ConfigurationManager cm, Corpus corpus) {
        this.cm = cm;
        this.corpus = corpus;

        buildUtteranceTree(corpus);

        buildWordTree(corpus);

        buildCharacterTree(corpus);

        buildPhonemeTree(corpus);
    }

    private void buildPhonemeTree(Corpus corpus) {
        DefaultMutableTreeNode wwn = new DefaultMutableTreeNode("Phonemes");

        for (String s : new TreeSet<String>( corpus.getPhonemeSequences())) {
            DefaultMutableTreeNode sn = new DefaultMutableTreeNode(s);
            for (Word w : new TreeSet<Word>(corpus.phonemeSequence2Words(s))) {
                DefaultMutableTreeNode wn = new DefaultMutableTreeNode( w );
                sn.add(wn);
            }
            wwn.add(sn);
        }

        phonemeTree = new JTree(wwn);
          phonemeTree.setFont(font);
        phonemeTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        phonemeTree.addTreeSelectionListener(this);

        phonemeScroll.setViewportView(phonemeTree);
    }

    private void buildCharacterTree(Corpus corpus) {

        DefaultMutableTreeNode wwn = new DefaultMutableTreeNode("Characters");

        for (String s : new TreeSet<String>( corpus.getCharacters()) ) {
            DefaultMutableTreeNode sn = new DefaultMutableTreeNode(hex2Unicode(s));
            for (Word w : new TreeSet<Word>( corpus.character2Words(s) )) {
                DefaultMutableTreeNode wn = new DefaultMutableTreeNode( w );
                sn.add(wn);
            }
            wwn.add(sn);
        }

        characterTree = new JTree(wwn);
        characterTree.setFont(font);
        characterTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        characterTree.addTreeSelectionListener(this);

        characterScroll.setViewportView(characterTree);
    }

    private void buildWordTree(Corpus corpus) {
        DefaultMutableTreeNode wwn = new DefaultMutableTreeNode("Words");

        for (String s : new TreeSet<String>(corpus.getSpellings())) {
            DefaultMutableTreeNode sn = new DefaultMutableTreeNode(hex2Unicode(s));
            for (Word w : corpus.getWords(s)) {
                DefaultMutableTreeNode wn = new DefaultMutableTreeNode( w );
                sn.add(wn);
            }
            wwn.add(sn);
        }

        wordTree = new JTree(wwn);
        wordTree.setFont(font);
        wordTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        wordTree.addTreeSelectionListener(this);

        wordScroll.setViewportView(wordTree);
    }

    private void buildUtteranceTree(Corpus corpus) {
        DefaultMutableTreeNode cn = new DefaultMutableTreeNode("Corpus and stuff");

        for (Utterance u : corpus.getUtterances()) {
            DefaultMutableTreeNode un = new DefaultMutableTreeNode(u.getRegionOfAudioData().getAudioDatabase().getPcm().getPcmFile() + " " + u.getBeginTime() + " " + u.getEndTime() + " " + u.getTranscript());
            for (Word w : u.getWords()) {
                DefaultMutableTreeNode wn = new DefaultMutableTreeNode( w );
                un.add(wn);
            }
            cn.add(un);
        }

        utteranceTree = new JTree(cn);
        utteranceTree.setFont(font);
        utteranceTree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        utteranceTree.addTreeSelectionListener(this);

        utteranceScroll.setViewportView(utteranceTree);
    }

    /**
     * Called whenever the value of the selection changes.
     *
     * @param e the event that characterizes the change.
     */
    public void valueChanged(TreeSelectionEvent e) {

        DefaultMutableTreeNode n = (DefaultMutableTreeNode) e.getPath().getLastPathComponent();
        if (n.getUserObject() instanceof Word) {
            Word word = (Word) n.getUserObject();

            System.out.println(word);


            final WordBrowser w = new WordBrowser(cm,word);
            JFrame f = new JFrame("WordBrowser");

            f.setContentPane(w.mainPane);
            f.setTitle(word.toString());
            f.pack();
            f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

            f.setVisible(true);
        }


    }

    public static String hex2Unicode(String hex) {
        if (hex.startsWith("<")) return hex;
        byte[] bytes = BatchNISTRecognizer.hex2Binary(hex);
        try {
            return new String(bytes, "GB2312");
        } catch (UnsupportedEncodingException e) {
            throw new Error(e);
        }
    }

    {
// GUI initializer generated by IntelliJ IDEA GUI Designer
// >>> IMPORTANT!! <<<
// DO NOT EDIT OR ADD ANY CODE HERE!
        $$$setupUI$$$();
    }

    /**
     * Method generated by IntelliJ IDEA GUI Designer
     * >>> IMPORTANT!! <<<
     * DO NOT edit this method OR call it in your code!
     */
    private void $$$setupUI$$$() {
        mainPanel = new JPanel();
        mainPanel.setLayout(new GridLayoutManager(1, 5, new Insets(0, 0, 0, 0), -1, -1));
        utteranceScroll = new JScrollPane();
        mainPanel.add(utteranceScroll, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null));
        utteranceScroll.setBorder(BorderFactory.createTitledBorder("utterances"));
        utteranceTree = new JTree();
        utteranceScroll.setViewportView(utteranceTree);
        wordScroll = new JScrollPane();
        mainPanel.add(wordScroll, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null));
        wordScroll.setBorder(BorderFactory.createTitledBorder("words"));
        wordTree = new JTree();
        wordScroll.setViewportView(wordTree);
        characterScroll = new JScrollPane();
        mainPanel.add(characterScroll, new GridConstraints(0, 2, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null));
        characterScroll.setBorder(BorderFactory.createTitledBorder("characters"));
        characterTree = new JTree();
        characterScroll.setViewportView(characterTree);
        phonemeScroll = new JScrollPane();
        mainPanel.add(phonemeScroll, new GridConstraints(0, 3, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null));
        phonemeScroll.setBorder(BorderFactory.createTitledBorder("phonemes"));
        phonemeTree = new JTree();
        phonemeScroll.setViewportView(phonemeTree);
        final JPanel panel1 = new JPanel();
        panel1.setLayout(new GridLayoutManager(2, 2, new Insets(0, 0, 0, 0), -1, -1));
        mainPanel.add(panel1, new GridConstraints(0, 4, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null));
        panel1.setBorder(BorderFactory.createTitledBorder("defaults"));
        final JPanel panel2 = new JPanel();
        panel2.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(panel2, new GridConstraints(1, 0, 1, 2, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_FIXED, GridConstraints.SIZEPOLICY_FIXED, null, null, null));
        final Spacer spacer1 = new Spacer();
        panel2.add(spacer1, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_VERTICAL, 1, GridConstraints.SIZEPOLICY_WANT_GROW, null, null, null));
        final JPanel panel3 = new JPanel();
        panel3.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(panel3, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null));
        panel3.setBorder(BorderFactory.createTitledBorder("bandwidth"));
        bandwidth = new JSlider();
        bandwidth.setOrientation(1);
        panel3.add(bandwidth, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null));
        final JPanel panel4 = new JPanel();
        panel4.setLayout(new GridLayoutManager(1, 1, new Insets(0, 0, 0, 0), -1, -1));
        panel1.add(panel4, new GridConstraints(0, 1, 1, 1, GridConstraints.ANCHOR_CENTER, GridConstraints.FILL_BOTH, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, GridConstraints.SIZEPOLICY_CAN_SHRINK | GridConstraints.SIZEPOLICY_CAN_GROW, null, null, null));
        panel4.setBorder(BorderFactory.createTitledBorder("offset"));
        offset = new JSlider();
        offset.setOrientation(1);
        offset.setPaintLabels(false);
        panel4.add(offset, new GridConstraints(0, 0, 1, 1, GridConstraints.ANCHOR_WEST, GridConstraints.FILL_HORIZONTAL, GridConstraints.SIZEPOLICY_WANT_GROW, GridConstraints.SIZEPOLICY_FIXED, null, null, null));
    }
}
