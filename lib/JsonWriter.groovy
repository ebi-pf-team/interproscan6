import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonException
import java.io.File
import java.io.FileWriter
import java.io.IOException

class JsonWriter {
    static void stream(String filePath, ObjectMapper mapper, Closure closure) {
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
            throw new JsonException("IO error writing file: $filePath", e)
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
