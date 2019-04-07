import std.stdio;
import std.range;
import std.algorithm;
import std.string;
import std.conv;
import std.file;
import std.path;
import std.parallelism;
import std.algorithm.comparison : equal;
import std.container.rbtree;
import std.traits;
import core.thread;
import std.process;
import std.datetime.stopwatch : benchmark;

struct FaiRecord {
    string header, lineTerm;
    ulong seqLen, lineLen, offset;
    @property ulong lineOffset() {
        return lineLen + lineTerm.length;
    }

    string toString() {
        return format("%s\t%s\t%s\t%s\t%s", header, seqLen, offset, lineLen, lineOffset);
    }
    unittest {
        auto rec = FaiRecord("chr2", "\n", 10, 50, 4);
        assert(rec.toString() == "chr2\t10\t4\t50\t51");
        rec.lineTerm = "\r\n";
        assert(rec.toString() == "chr2\t10\t4\t50\t52");
    }

    this(string str) {
        auto res = str.split("\t");
        header ~= res[0];
        seqLen = to!ulong(res[1]);
        offset = to!ulong(res[2]);
        lineLen = to!ulong(res[3]);
        lineTerm = (to!ulong(res[4])-lineLen) == 1 ? "\n" : "\r\n";
    }
    unittest {
        auto s = "chr2\t10\t4\t50\t51";
        assert(FaiRecord(s).toString() == s);
    }

    this(string header, string lineTerm, ulong seqLen, ulong lineLen, ulong offset) {
        this.header = header;
        this.seqLen = seqLen;
        this.offset = offset;
        this.lineLen = lineLen;
        this.lineTerm = lineTerm;
    }
    unittest {
        assert(FaiRecord("chr2", "\n", 10, 50, 4).toString() == "chr2\t10\t4\t50\t51");
    }
}


auto buildFai() {

    File f = File("/home/ubuntu/Documents/Dexp/TEST/fastaExmp.fna", "r"); //f.seek(0); file may be gone, exception?
    string lineTerm = f.byLine(KeepTerminator.yes).take(1).front.endsWith("\r\n") ? "\r\n" : "\n";

    
    ulong lineLen = f.byLine(KeepTerminator.yes).take(1).to!string.length - 4; //to!string adds brackets like ["%data%"]
    //f.byLine(KeepTerminator.yes).take(1).writeln();
    
    //Maybe try block?
    string[] contents;
    try{
        contents = readText("/home/ubuntu/Documents/Dexp/TEST/fastaExmp.fna").splitLines(KeepTerminator.no);
    }
    catch(FileException e){
        return parseFastaSingleThread();
    }
    FaiRecord[] records; 
    auto readHeaders = redBlackTree!string;

    ulong step = f.size()/taskPool.size();
    Task!(parseFasta, Parameters!parseFasta)*[] tasks;//(parseFasta, Parameters!parseFasta)*[] tasks;
    // //writeln("File Size - ", f.size(), "\ntaskPool Size - ", taskPool.size(), "\nlineLen - ", lineLen, "\nstep - ", step );

    //writeln(taskPool.size());
    for(int i = 0; i < taskPool.size(); i++){ 
        tasks ~= task!parseFasta(&contents, &readHeaders, i*step/lineLen, lineTerm);
        tasks[i].executeInNewThread();
    }

    foreach(t; tasks){
        records ~= t.yieldForce();
    }

    return records;

}

auto parseFasta(string[]* contents, RedBlackTree!(string)* readHeaders, ulong firstLineNumber, string lineTerm){

    FaiRecord[] records; 
    ulong offset;

    
    foreach(line; (*contents)[firstLineNumber..$]) {
        writeln('(',line,'|',(*contents).back,')');
        offset += line.length + lineTerm.length;
        if ( line.startsWith(">") ) {
			string header = line.split(" ").front[1..$];
            
            if(header in (*readHeaders)){
                return records;
            }
            (*readHeaders).insert(header);

            records~=FaiRecord();
            records[$-1].lineTerm = lineTerm;
            records[$-1].header ~= header;
            records[$-1].offset = offset;
        } else {
            if(records.length == 0) continue;
            //writeln("Records Length - ", records.length, "\n Line - ", line,  "\n Thread ID - ", thisThreadID(),']');
            records[$-1].toString();
            if ( records[$-1].lineLen == 0 ) {
				records[$-1].lineLen = line.length;
            }
            records[$-1].seqLen += line.length;
        }
    }

    return records;

}

auto parseFastaSingleThread(){
    File f = File("/home/ubuntu/Documents/Dexp/TEST/fastaExmp.fna", "r");
    FaiRecord[] records; 
    string lineTerm = f.byLine(KeepTerminator.yes).take(1).front.endsWith("\r\n") ? "\r\n" : "\n";
    f.seek(0);
    ulong offset;
    foreach(line; f.byLine(KeepTerminator.no, lineTerm)) {
        offset+= line.length + lineTerm.length;
        if ( line.startsWith(">") ) {
            records~=FaiRecord();
            records[$-1].lineTerm = lineTerm;
            records[$-1].header ~= line.split(" ").front[1..$];
            records[$-1].offset = offset;
        } else {
            if ( records[$-1].lineLen == 0 ) {
                records[$-1].lineLen = line.length;
            }
            records[$-1].seqLen += line.length;
        }
    }

    return records;
}



int main(){
    //auto records = parseFastaSingleThread("/home/ubuntu/Documents/Dexp/TEST/fasta");//buildFai("/home/ubuntu/Documents/Dexp/TEST/fasta");
    //foreach(t; records) writeln(t.toString());
    auto r = benchmark!(buildFai, parseFastaSingleThread)(1);
    writeln("Parallel - ",r[0],"\nSingle thread - ",r[1]);
    return 0;
}