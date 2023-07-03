<script>
function buildQuiz(myq, qc){
  // variable to store the HTML output
  const output = [];

  // for each question...
  myq.forEach(
    (currentQuestion, questionNumber) => {

      // variable to store the list of possible answers
      const answers = [];

      // and for each available answer...
      for(letter in currentQuestion.answers){

        // ...add an HTML radio button
        answers.push(
          `<label>
            <input type="radio" name="question${questionNumber}" value="${letter}">
            ${letter} :
            ${currentQuestion.answers[letter]}
          </label><br/>`
        );
      }

      // add this question and its answers to the output
      output.push(
        `<div class="question"> ${currentQuestion.question} </div>
        <div class="answers"> ${answers.join('')} </div><br/>`
      );
    }
  );

  // finally combine our output list into one string of HTML and put it on the page
  qc.innerHTML = output.join('');
}

function showResults(myq, qc, rc){

  // gather answer containers from our quiz
  const answerContainers = qc.querySelectorAll('.answers');

  // keep track of user's answers
  let numCorrect = 0;

  // for each question...
  myq.forEach( (currentQuestion, questionNumber) => {

    // find selected answer
    const answerContainer = answerContainers[questionNumber];
    const selector = `input[name=question${questionNumber}]:checked`;
    const userAnswer = (answerContainer.querySelector(selector) || {}).value;

    // if answer is correct
    if(userAnswer === currentQuestion.correctAnswer){
      // add to the number of correct answers
      numCorrect++;

      // color the answers green
      answerContainers[questionNumber].style.color = 'lightgreen';
    }
    // if answer is wrong or blank
    else{
      // color the answers red
      answerContainers[questionNumber].style.color = 'red';
    }
  });

  // show number of correct answers out of total
  rc.innerHTML = `${numCorrect} out of ${myq.length}`;
}
</script>


The dataset used in this workshop is a subset of a much larger dataset from a [recent study](https://doi.org/10.1038/s41588-022-01088-x) that generated single nuclei transcriptome and chromatin accessibility profiles from colorectal tissue samples. The authors isolated 1000 to 10000 nuclei each from 81 samples of three types: 48 polyp samples, 27 normal tissue samples, and 6 colorectal cancer (CRC) samples from patients with or without germline APC mutations. They observed a continuum of cell state and composition changes from normal tissue, to polyps, to cancer.

For the purposes of this workshop, we will use one sample from each condition (CRC: A001-C-007, polyp: A001-C-104, and normal: B001-A-301).

Source:

Becker, W. R.; Nevins, S. A.; Chen, D. C.; Chiu, R.; Horning, A. M.; Guha, T. K.; Laquindanum, R.; Mills, M.; Chaib, H.; Ladabaum, U.; Longacre, T.; Shen, J.; Esplin, E. D.; Kundaje, A.; Ford, J. M.; Curtis, C.; Snyder, M. P.; Greenleaf, W. J. Single-Cell Analyses Define a Continuum of Cell State and Composition Changes in the Malignant Transformation of Polyps to Colorectal Cancer. Nat. Genet. 2022. https://doi.org/10.1038/s41588-022-01088-x.

# Data Setup

Let's set up a project directory for the analysis, and talk a bit about project philosophy.

**1\.** First, create a directory for your user and the example project in the workshop directory:

```bash
cd
mkdir -p /share/workshop/scRNA_workshop/$USER/scrnaseq_example
```

---

**2a\.** Next, go into that directory, create a raw data directory (we are going to call this 00-RawData) and cd into that directory. Let's then create symbolic links to the fastq files that contains the raw read data.

```bash
cd /share/workshop/scRNA_workshop/$USER/scrnaseq_example
mkdir 00-RawData
cd 00-RawData/
ln -s /share/workshop/scRNA_workshop/Data/*.fastq.gz .
```

This directory now contains the reads for each sample.

**2b\.** Let's create a sample sheet for the project, and store sample names in a file called samples.txt

```bash
cd /share/workshop/scRNA_workshop/$USER/scrnaseq_example/00-RawData
ls *_R1_* |cut -d'_' -f1 - > ../samples.txt
cat ../samples.txt
```

# Data Exploration
---
**3\.** Now, take a look at the raw data directory.

```bash
ls /share/workshop/scRNA_workshop/$USER/scrnaseq_example/00-RawData
```
---

**4\.** View the contents of the files using the 'less' command, when gzipped used 'zless' (which is just the 'less' command for gzipped files, q to exit):

Read 1

```bash
cd /share/workshop/scRNA_workshop/$USER/scrnaseq_example/00-RawData
zless A001-C-007_S4_R1_001.fastq.gz
```

and Read 2

```bash
zless A001-C-007_S4_R2_001.fastq.gz
```

A detailed explanation of FASTQ file can be found [here](filetypes.md). Please read on the description and make sure you can identify which lines correspond to a single read and which lines are the header, sequence, and quality values. Press 'q' to exit this screen.

Let's figure out the number of reads in this file. A simple way to do that is to count the number of lines and divide by 4 (because the record of each read uses 4 lines). In order to do this use cat to output the uncompressed file and pipe that to "wc" to count the number of lines:

```bash
zcat A001-C-007_S4_R1_001.fastq.gz | wc -l
```

Divide this number by 4 and you have the number of reads in this file. One more thing to try is to figure out the length of the reads without counting each nucleotide. First get the first 4 lines of the file (i.e. the first record):

```bash
zcat A001-C-007_S4_R1_001.fastq.gz  | head -4
```

Note the header lines (1st and 3rd line) and sequence and quality lines (2nd and 4th) in each 4-line fastq block. You can isolate the sequence line:

```bash
zcat A001-C-007_S4_R1_001.fastq.gz | head -2 | tail -1
```

Then, copy and paste the DNA sequence line into the following command (replace [sequence] with the line):

```bash
echo -n [sequence] | wc -c
```

This will give you the length of the read.

Also can do the bash one liner:

```bash
echo -n $(zcat A001-C-007_S4_R1_001.fastq.gz  | head -2 | tail -1) | wc -c
```

See if you can figure out how this command works.


## Quiz

<div id="quiz1" class="quiz"></div>
<button id="submit1">Submit Quiz</button>
<div id="results1" class="output"></div>
<script>
quizContainer1 = document.getElementById('quiz1');
resultsContainer1 = document.getElementById('results1');
submitButton1 = document.getElementById('submit1');

myQuestions1 = [
  {
    question: "How many reads are in the file?",
    answers: {
      a: "200 Thousand",
      b: "500 Thousand",
      c: "1 Million",
      d: "2 Million"
    },
    correctAnswer: "c"
  },
  {
    question: "What is the length of Read 1 and Read 2?",
    answers: {
      a: "90 and 90",
      b: "50 and 100",
      c: "28 and 91",
      d: "75 and 25"
    },
    correctAnswer: "c"
  }
];

buildQuiz(myQuestions1, quizContainer1);
submitButton1.addEventListener('click', function() {showResults(myQuestions1, quizContainer1, resultsContainer1);});
</script>
