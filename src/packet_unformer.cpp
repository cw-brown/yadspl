#include "packet_unformer.h"
#include "crcinterface.h"
#include <bit>

PacketUnformer::PacketUnformer(
    std::vector<bool> * in_data_buffer,
    std::vector<uint8_t> *  out_data_buffer) : 
    inDataBuffer(in_data_buffer),
    outDataBuffer(out_data_buffer) {

    // Make crc generator with polynomial for ISO 3309 in reversed format
    crcGenny = crcutil_interface::CRC::Create(0xEDB88320, 0, 32, true, 0, 0, 0, true, NULL);

}

uint8_t PacketUnformer::formPacket() {

    // Temp data storage
    uint8_t tempFlag;
    uint8_t tempLen;
    bool noFlag;

    // Setup
    noFlag = true;
    inputPos = inDataBuffer->begin();

    // Find flag
    while(noFlag) {

        tempFlag = 0;

        // First read potential flag field, checking buffer has enough data
        for(int i = 0; i < 8; i++) {

            if(inputPos + i == inDataBuffer->cend()) {

                return 0;

            }

            tempFlag += (*(inputPos + i) << (7 - i));

        }

        // Check if candidate flag good
        noFlag = (tempFlag != 0b01111110);

        // Move up pointer
        inputPos++;

        // If buffer empty, return
        if(inputPos == inDataBuffer->cend()) {

            return 0;

        }

    }

    // Check for second flag


    // Read potential len field: inputPos one bit after flag, need to check 4 forward for length, checking buffer has enough data
    for(int i = 0; i < 4; i++) {

        if(inputPos + i == inDataBuffer->cend()) {

            return 0;
            
        }

        tempLen += (*((inputPos+4) + i) << (3 - i));

    }

    // Second flag should be 8 + 64 + 8 + 8*length + 32 bits ahead of inputPos

    // Check that enough data is present in buffer
    for(int i = 0; i < 8 + 64 + 8 + 8*tempLen + 32; i++) {

        if(inputPos + i == inDataBuffer->cend()) {

            return 0;
            
        }

    }

    // Read potential flag field
    for(int i = 0; i < 8; i++) {

        tempFlag += (*((inputPos + 8 + 64 + 8 + 8*tempLen + 32) + i) << (7 - i));

    }

    if(tempFlag != 0b01111110) {

        return 1;

    }

}