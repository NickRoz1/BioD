
import bio.std.hts.bam.reader;
import bio.std.hts.bam.baseinfo;

import std.stdio;
import std.range : take, drop;
import std.algorithm : find;

void main() {

  auto bam = new BamReader("../test/data/b7_295_chunk.bam");

  // get read group information by name
  auto rg = bam.header.read_groups["9IKNG"];

  auto read = find!(r => r.name == "9IKNG:00592:01791")(bam.reads).front;

  // fetch information about flow calls from FZ & ZF tags
  // and also reference base from MD tag
  auto bases = basesWith!("FZ", "MD")(read, arg!"flowOrder"(rg.flow_order),
      arg!"keySequence"(rg.key_sequence));

  // end of read contains a few indel errors
  foreach (baseinfo; bases.drop(350).take(32)) {
    writefln("%s\t%s\tflow: %3d\tintensity: %.2f\t\tref. pos.: %6d\tCIGAR op.: %s", 
        baseinfo.reference_base,        // from MD tag
        baseinfo.base,
        baseinfo.flow_call.flow_index,  // from FZ tag
        baseinfo.flow_call.intensity,   // also from FZ tag
        baseinfo.position, 
        baseinfo.cigar_operation);
    // notice that because the read is on reverse strand, 
    // reference position decreases during the iteration
  }
}
