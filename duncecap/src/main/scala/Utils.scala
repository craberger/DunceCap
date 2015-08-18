package DunceCap

import scala.collection.mutable
import sys.process._
import scala.io._
import java.nio.file.{Paths, Files}
import java.io.{FileWriter, File, BufferedWriter}

class CodeStringBuilder {
  val buffer = new StringBuilder
  def println(str : String): Unit = {
    buffer.append(str + "\n")
  }
  def print(str : String): Unit = {
    buffer.append(str)
  }

  override def toString = {
    buffer.toString
  }
}

object Utils{
  def loadEnvironmentFromJSON(config:Map[String,Any], create_db:Boolean):Unit = {
    //Pull top level vars from JSON file
    val db_folder = config("database").asInstanceOf[String]
    Environment.layout = config("layout").asInstanceOf[String]
    Environment.algorithm = config("algorithm").asInstanceOf[String]
    Environment.numNUMA = config("numNUMA").asInstanceOf[Double].toInt
    Environment.numThreads = config("numThreads").asInstanceOf[Double].toInt
    Environment.dbPath = db_folder

    if(create_db){
      //setup folders for the database on disk
      s"""mkdir -p ${db_folder} ${db_folder}/relations ${db_folder}/encodings""" !
    }
    
    val relations = config("relations").asInstanceOf[List[Map[String,Any]]]
    relations.foreach(r => {
      val name = r("name").asInstanceOf[String]
      val source = r("source").asInstanceOf[String]
      val attributes = r("attributes").asInstanceOf[List[Map[String,String]]]

      attributes.foreach(a => { 
        Environment.addEncoding(new Encoding(a("encoding"),a("type")))
        if(create_db)
          s"""mkdir -p ${db_folder}/encodings/${a("encoding")} """ !
      })
      val attribute_names = attributes.map(a => a("name"))
      val attribute_types = attributes.map(a => a("type"))
      val attribute_encodings = attributes.map(a => a("encoding"))

      val ordering = r("ordering").asInstanceOf[List[String]]
      val orderings = if(ordering.length == 1 && ordering(0) == "all") attribute_names.permutations.toList else List(ordering)
      orderings.foreach(o => {      
        val o_name = o.mkString("")
        val rel_name = s"""${name}_${o_name}""" 
        val relationIn = new Relation(rel_name,attribute_names,attribute_types,attribute_encodings)
        Environment.addRelation(name,relationIn)
        if(create_db){
          Environment.addASTNode(ASTLoadRelation(source,relationIn))
          s"""mkdir ${db_folder}/relations/${rel_name}""" !
        }
      })
    })
  }

  def compileAndRun(codeStringBuilder:CodeStringBuilder):Unit = {
    val cppDir = "../emptyheaded/generated"
    val dbName = Environment.dbPath.split("/").toList.last
    val cppFilepath = s"""${cppDir}/load_${dbName}.cpp"""
    if (!Files.exists(Paths.get(cppDir))) {
      println(s"""Making directory ${cppDir}""")
      s"""mkdir ${cppDir}""" !
    } else {
      println(s"""Found directory ${cppDir}, continuing...""")
    }

    val file = new File(cppFilepath)
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(codeStringBuilder.toString)
    bw.close()
    s"""clang-format -style=llvm -i ${file}""" !

    sys.process.Process(Seq("rm", "-rf",s"""bin/load_${dbName}"""), new File("../emptyheaded")).!
    val result = sys.process.Process(Seq("make",s"""NUM_THREADS=${Environment.numThreads}""", s"""bin/load_${dbName}"""), new File("../emptyheaded")).!

    if (result != 0) {
      println("FAILURE: Compilation errors.")
      System.exit(1)
    }

    s"""../emptyheaded/bin/load_${dbName}""" !
  }
}