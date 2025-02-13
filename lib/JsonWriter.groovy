import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonException
import java.io.File
import java.io.FileWriter
import java.io.IOException

class JsonWriter {
    static void streamJson(String filePath, ObjectMapper mapper, Closure closure) {
        /* Write out json objects, e.g. at the end of XREFS
        To use:
        ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
        JsonWriter.streamMap(outputFilePath.toString(), jacksonMapper) { JsonGenerator jsonGenerator ->
            jsonGenerator.writeStringField("exampleKey", "exampleValue")
        }
        */
        FileWriter fileWriter = null
        JsonGenerator generator = null
        try {
            JsonFactory factory = mapper.getFactory()
            fileWriter = new FileWriter(new File(filePath))
            generator = factory.createGenerator(fileWriter)
            generator.writeStartObject()

            closure.call(generator)  // Call the closure to write key-value pairs

            generator.writeEndObject()
        } catch (IOException e) {
            throw new JsonException("IO error writing file: $filePath -- $e\n${e.printStackTrace()}", e)
        } catch (Exception e) {
            throw new Exception("Error occured when writing Json file $filePath -- $e\n${e.printStackTrace()}", e)
        } finally {
            if (generator != null) {
                generator.close()
            }
            if (fileWriter != null) {
                fileWriter.close()
            }
        }
    }

    static void streamArray(String filePath, ObjectMapper mapper, Closure closure) {
        // For writing out Json Arrays, e.g. at the end of Aggergage_All_Matches
        FileWriter fileWriter = null
        JsonGenerator generator = null
        try {
            JsonFactory factory = mapper.getFactory()
            fileWriter = new FileWriter(new File(filePath))
            generator = factory.createGenerator(fileWriter)
            generator.writeStartArray()
            closure.call(generator)
            generator.writeEndArray()
        } catch (IOException e) {
            throw new JsonException("IO error writing file: $filePath -- $e\n${e.printStackTrace()}", e)
        } finally {
            if (generator != null) {
                generator.close()
            }
            if (fileWriter != null) {
                fileWriter.close()
            }
        }
    }

    static void writeMap(JsonGenerator generator, ObjectMapper mapper, Object data) throws IOException {
        /* Used for writing out a Map object in the middle of a Json writer stream.
        For example during XREFS it is used to write out each Match to the XREFs output JSON file.
        To use:
        ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
        JsonWriter.streamMap(outputFilePath.toString(), jacksonMapper) { JsonGenerator jsonGenerator ->
            JsonWriter.writeMap(jsonGenerator, jacksonMapper, Match.asMap(match))
         }
         */
        if (data instanceof Map) {
            data.each { key, value ->
                generator.writeFieldName(key.toString())
                if (value instanceof Map) {
                    generator.writeStartObject()
                    value.each { k, v ->
                        generator.writeFieldName(k.toString())
                        mapper.writeValue(generator, v)
                    }
                    generator.writeEndObject()
                } else {
                    mapper.writeValue(generator, value)
                }
            }
        } else {
            mapper.writeValue(generator, data)
        }
    }

    static void writeMaptoFile(String filePath, ObjectMapper mapper, Map data) {
        /* Method to write (i.e. dump) a map or object (like seqMatchesAggreg) to a file as a whole, i.e.
        when not streaming.
        */
        FileWriter fileWriter = null
        JsonGenerator generator = null
        try {
            JsonFactory factory = mapper.getFactory()
            fileWriter = new FileWriter(new File(filePath))
            generator = factory.createGenerator(fileWriter)
            generator.writeStartObject()

            data.each { key, value ->
                generator.writeFieldName(key)
                mapper.writeValue(generator, value)
            }

            generator.writeEndObject()
        } catch (IOException e) {
            throw new JsonException("IO error writing file: $filePath -- $e\n${e.printStackTrace()}", e)
        } finally {
            if (generator != null) {
                generator.close()
            }
            if (fileWriter != null) {
                fileWriter.close()
            }
        }
    }
}
