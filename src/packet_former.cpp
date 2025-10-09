#include "packer_former.h"
#include "crcinterface.h"

PacketFormer::PacketFormer(
        uint32_t sender_address,
        uint32_t reciever_address,
        uint8_t data_length) :
        //inputBuffer(input_buffer),
        //outputBuffer(output_buffer),
        senderAddress(sender_address),
        recieverAddress(reciever_address),
        dataLength(data_length) {

    // Make crc generator with polynomial for ISO 3309 in reversed format
    crcGenny = crcutil_interface::CRC::Create(0xEDB88320, 0, 32, true, 0, 0, 0, true, NULL);
    tempPacket.senderAddress = senderAddress;
    tempPacket.recieverAddress = recieverAddress;

}



void PacketFormer::setDataLength(uint8_t data_length) {

    dataLength = data_length;

}

Packet PacketFormer::formDataPacket(std::vector<uint8_t> * data, uint8_t sequence_number) {

    crcutil_interface::UINT64 tempCRC;

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    
    tempPacket.controlCode = 0b0000;

    tempPacket.dataLength = dataLength;

    tempPacket.sequenceNumber = sequence_number;

    // Form data section
    for(int i = dataLength - 1; i >= 0; i--) {

        (tempPacket.data).push_back(*data->end());
        data->pop_back();

    }

    // Generate CRC and store to ERC field
    crcGenny->Compute(tempPacket.data.data(), tempPacket.dataLength, &tempCRC, NULL);
    tempPacket.erc = (uint32_t) tempCRC;

    return tempPacket;

}

Packet PacketFormer::formRetransmitPacket(std::vector<uint8_t> * data, uint8_t sequence_num) {
    
    // Reset packet contents
    (tempPacket.data).clear();
    
    // Form initial packet with class method
    tempPacket = formDataPacket(data, sequence_num);

    // Set packet field(s)

    tempPacket.controlCode = 0b0001;

    return tempPacket;

}

Packet PacketFormer::formRetransmitPacket(Packet data_packet) {

    // Set packet field(s)
    data_packet.controlCode = 0b0001;

    return data_packet;

}

Packet PacketFormer::formBusyStartPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0010;

    tempPacket.dataLength = 0;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

}

Packet PacketFormer::formBusyEndPacket() {

    // Reset packet contents
    (tempPacket.data).clear();

    // Set packet field(s)
    tempPacket.controlCode = 0b0011;

    tempPacket.dataLength = 0;

    tempPacket.sequenceNumber = 0;

    tempPacket.erc = 0;

}