import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonException
import java.io.File
import java.io.FileWriter
import java.io.IOException

class JsonWriter {
    static void streamJson(String filePath, ObjectMapper mapper, Closure closure) {
        /* Write out json objects while streaming the file, e.g. when writing the final JSON file
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
            throw new JsonException("IO error writing file: $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
        } catch (Exception e) {
            throw new Exception("Error occured when writing Json file $filePath -- $e\n${e.printStackTrace()}\n${e.getCause()}", e)
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
