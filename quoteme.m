function quote = quoteme()
% QUOTEME   an inspirational quote for your code!
%   quote = QUOTEME delivers an beautiful (or snarky) quote about science
%   from history's best. 
%
%   Last updated 10/14/2016

Quotes{1} = 'If your experiment needs statistics, you ought to have done a better experiment.\nErnest Rutherford, 1871 to 1937';
Quotes{2} = 'An experiment is a question which science poses to Nature, and a measurement is the recording of Nature’s answer.\nMax Planck, 1858 to 1947';
Quotes{3} = 'Falsity in intellectual action is intellectual immorality.\nThomas Chrowder Chamberlin, 1843 to 1928';
Quotes{4} = 'It does not help that some politicians and journalists assume the public is interested only in those aspects of science that promise immediate practical applications to technology or medicine.\nSteven Weinberg, 1933 to present';
Quotes{5} = 'Valid criticism does you a favor.\nCarl Sagan, 1934 to 1996';
Quotes{6} = 'Pierre Curie voluntarily exposed his arm to the action of radium for several hours. This resulted in damage resembling a burn that developed progressively and required several months to heal. Henri Becquerel had by accident a similar burn as a result of carrying in his vest pocket a glass tube containing radium salt. He came to tell us of this evil effect of radium, exclaiming in a manner at once delighted and annoyed: “I love it, but I owe it a grudge.”\nMarie Curie, 1867 to 1934';
Quotes{7} = 'I believe there are 15 747 724 136 275 002 577 605 653 961 181 555 468 044 717 914 527 116 709 366 231 425 076 185 631 031 296 protons in the universe and the same number of electrons.\nSir Arthur Eddington, 1882 to 1944';
Quotes{8} = 'God created two acts of folly. First, He created the Universe in a Big Bang. Second, He was negligent enough to leave behind evidence for this act, in the form of microwave radiation.\nPaul Erdős, 1913 to 1996';
Quotes{9} = 'Nothing in life is to be feared. It is only to be understood.\nMarie Curie (1867-1934)';
Quotes{10} = 'Science and everyday life cannot and should not be separated.\nRosalind Franklin (1920-1958)';
Quotes{11} = 'The Sun, with all the planets revolving around it, and depending on it, can still ripen a bunch of grapes as though it had nothing else in the Universe to do.\nGalileo Galilei, 1564-1642';
Quotes{12} = 'A scientist in his laboratory is not a mere technician; he is also a child confronting natural phenomena that impress him as though they were fairy tales.\nMarie Curie, 1867-1932';
Quotes{13} = 'A ship in port is safe, but that is not what ships are for. Sail out to sea and do new things.\nGrace Hopper, 1906-1992';
Quotes{14} = 'There are an awful lot of scientists today who believe that before very long we shall have unraveled all the secrets of the universe. There will be no puzzles anymore. To me it’d be really, really tragic because I think one of the most exciting things is this feeling of mystery, feeling of awe, the feeling of looking at a little live thing and being amazed by it and how it’s emerged through these hundreds of years of evolution, and there it is, and it is perfect and why.\nJane Goodall, 1934-';
Quotes{15} = 'So my antagonist said, "Is it impossible that there are flying saucers? Can you prove that it is impossible?" "No", I said, "I cannot prove it is impossible. It is just very unlikely". At that he said, "You are very unscientific. If you cannot prove it impossible then how can you say that it is unlikely?" But that is the way that is scientific. It is scientific only to say what is more likely and what less likely, and not to be proving all the time the possible and impossible.\nRichard Feynman, 1918-1988';
Quotes{16} = 'What we know here is very little, but what we are ignorant of is immense.\nPierre Laplace, 1749-1829';
Quotes{17} = 'I like understanding things and explaining them, and sometimes when you’re trying to understand something, you see something new, and they call that research.\nDavid Blackwell, 1919-2010';
Quotes{18} = 'You look at science (or at least talk of it) as some sort of demoralising invention of man, something apart from real life, and which must be cautiously guarded and kept separate from everyday existence. But science and everyday life cannot and should not be separated.\nRosalind Franklin 1928-1958'; 


winner = 'CONGRATULATIONS!!! YOU"RE A WINNER!!! YOU ARE THE 100th PERSON TO RUN THIS CODE!!! Click here to claim your FREE prize - One completely bug-free code run sponsored by The Universe(TM)';

if randi(100) == 100;
    quote = wraptext(winner);
else
    quote = wraptext(sprintf(Quotes{randi(length(Quotes))}));
end
